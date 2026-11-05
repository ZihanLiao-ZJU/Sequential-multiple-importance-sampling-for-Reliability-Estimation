classdef SeMIS_low_dimension
    % SeMIS with Nataf transformation for structural reliability problems
    % Low-dimensional variant: Kriging-based hybrid truncation (Eqn. 25).
    % ISD (Eqn. 21): q_i(u) = phi_n(u; sigma^2 I) * L_i(u) / z_i
    %   rPDF pi_i(u) = phi_n(u; sigma^2 I)                (Eqn. 12)
    %   ISF  L_i(u)  = I[g(u)<l_i] * prod_m I[g_hat_m(u)<=l_hat_m]  (Eqn. 20)
    properties
        TF
        Nfun
        Ndim
        Ncal
        X
        Y
        Nsam = 1000;
    end

    properties
        G
    end

    properties(Dependent)
        FlgCvg
    end

    methods
        function obj = SeMIS_low_dimension(PF)
            % Initialize the SeMIS_low_dimension object
            obj.TF   = PF;
            obj.Nfun = PF.Nfun + 3;
            obj.Ndim = PF.Ndim;

            % Approximate optimal variance-inflation factor (Eqn. 19):
            %   sigma_tilde* = sqrt(1 + beta * lambda(beta) / n)
            % with beta = 2.33 (target P_f ~ 1e-2) and lambda(beta) the inverse Mills ratio.
            beta      = 2.33;
            lambda_beta = normpdf(beta) ./ (1 - normcdf(beta));
            sigma     = sqrt((obj.Ndim + beta * lambda_beta) ./ obj.Ndim);

            obj.G.p      = 0.1;      % intermediate level probability p (Section 4.1.2)
            obj.G.Ite    = 0;        % current iteration counter
            obj.G.l      = inf;      % exact PF threshold l_i (Eqn. 20)
            obj.G.g_hat  = {};       % Kriging surrogates g_hat_m (Eqn. 25)
            obj.G.l_hat  = [];       % surrogate thresholds l_hat_m (Eqn. 27)
            obj.G.kappa  = 2.0;      % conservatism factor kappa (Eqn. 25)

            % rPDF parameters (Eqn. 12): pi_i(u) = phi_n(u; sigma^2 I)
            obj.G.L_pi  = eye(obj.Ndim) * sigma;    % Cholesky factor (sigma*I)
            obj.G.Sigma = obj.G.L_pi * obj.G.L_pi.';% covariance sigma^2 * I
            obj.G.mu_pi = zeros(obj.Ndim,1);        % mean vector (0)

            obj.Ncal = 0;
        end

        function obj = EvlY(obj,u)
            % Evaluate the function values associated with the input samples
            u = u(1:end,:);
            obj.X = u;
            n_sam = size(u,2);

            % Cheap component used in the standard normal space: log phi_n(u)
            log_phi = logGauss(u);

            % Screening condition defined by Kriging in the U-space:
            %   g_hat_m(u) = mu_hat(u) - kappa * sd_hat(u)  (Eqn. 25)
            if obj.G.Ite > 0
                ind_scr = true(1,n_sam);
                for i = 1:obj.G.Ite
                    [mu_hat,sd_hat] = predict(obj.G.g_hat{i},u.');
                    g_hat = (mu_hat - obj.G.kappa .* sd_hat).';
                    ind_scr = ind_scr & (g_hat <= obj.G.l_hat(i));
                end
            else
                ind_scr = true(1,n_sam);
            end

            % Default setting: no PF evaluation is performed
            g_val = inf(1,n_sam);
            n_cal = 0;

            if any(ind_scr)
                x       = obj.TF.Dist.U2X(u(:,ind_scr));
                g_val(ind_scr) = obj.TF.EvlLSF(x);
                n_cal   = sum(ind_scr);
            end

            obj.Y    = [log_phi; g_val];   % [log phi_n(u); g(u)]
            obj.Ncal = n_cal;
        end

        function log_L = EvlLKF(obj)
            % Evaluate the intermediate log-likelihood = log ISF L_i(u) (Eqn. 20):
            %   log L_i(u) = log I[g(u)<l_i] + sum_m log I[g_hat_m(u)<=l_hat_m]
            g_val = obj.Y(2,:);
            l     = obj.G.l;
            n_sam = size(g_val,2);

            if isinf(l)
                log_L = zeros(1,n_sam);
                return;
            end

            if obj.G.Ite > 0
                ind_scr = true(1,n_sam);
                for i = 1:obj.G.Ite
                    [mu_hat,sd_hat] = predict(obj.G.g_hat{i},obj.X.');
                    g_hat = (mu_hat - obj.G.kappa .* sd_hat).';
                    ind_scr = ind_scr & (g_hat <= obj.G.l_hat(i));
                end
            else
                ind_scr = true(1,n_sam);
            end

            log_L = log(double((g_val < l) & ind_scr));
        end

        function log_pi = EvlPDF(obj)
            % Evaluate the intermediate log-prior density = log rPDF pi_i(u) (Eqn. 12)
            log_pi = logGauss(obj.X,obj.G.mu_pi,obj.G.Sigma);
        end

        function [y,g,obj,N_cal] = UpdObj(obj,u,y,p)
            % Update the intermediate distribution parameters
            %   l_i       : exact PF threshold = p-quantile of g (Eqn. 24)
            %   g_hat_{i+1}: Kriging surrogate trained on evaluated samples (Eqn. 25)
            %   l_hat     : conservative surrogate threshold = max over seeds (Eqn. 27)
            g_val = y(4,:);
            u = u(1:end,:);
            g = obj.G;
            if nargin >= 4, g.p = p; end

            l        = max(quantile(g_val,g.p),0);     % threshold l_i (Eqn. 24)
            ind_seed = g_val < l;                      % seed set S_{i+1} (Eqn. 24)

            % Fit kriging using all truly evaluated samples
            ind_fit = isfinite(g_val);
            u_fit   = u(:,ind_fit);
            g_fit   = g_val(ind_fit);

            mdl = fitrgp( ...
                u_fit.',g_fit.', ...
                'KernelFunction','ardsquaredexponential', ...
                'BasisFunction','constant', ...
                'Standardize',true, ...
                'FitMethod','exact', ...
                'PredictMethod','exact');

            % Conservative kriging score (Eqn. 25)
            [mu_hat,sd_hat] = predict(mdl,u.');
            g_hat = (mu_hat - g.kappa .* sd_hat).';

            % Threshold chosen so that all current seeds are retained (Eqn. 27)
            l_hat_tmp = max(g_hat(ind_seed));

            g.g_hat{end+1} = mdl;
            g.l_hat        = cat(1,g.l_hat,l_hat_tmp);
            g.l            = l;
            g.Ite = g.Ite + 1;

            obj.G   = g;
            N_cal   = 0;
        end

        function FlgCvg = get.FlgCvg(obj)
            % Check whether the convergence criterion is satisfied (Section 4.1.2)
            FlgCvg = obj.G.l <= 0;
        end
    end
    methods
        function defaults=run_defaults(obj)
            defaults=struct('NumSam',1000,'ConPrb',0.1,'MaxIte',10000,'Verbose',true);
            if isprop(obj,'Nsam'),defaults.NumSam=obj.Nsam;end
            if isfield(obj.G,'p'),defaults.ConPrb=obj.G.p;end
        end
        function c=start_run(~,owner,sampler)
            c=SeMCS.start_context(owner,sampler,owner.MaxIte-1,'replace');
        end
        function [c,r]=sample_level(~,c,level)
            [c,r]=SeMCS.sample_context(c,level);
        end
        function c=advance_level(~,c,r,level)
            if level+1<=numel(c.owner.G)
                c.model.G=SeMCS.prescribed_state(c.model.G,c.owner.G{level+1},false,'replace');extra=0;
            else
                [c.y{level},~,c.model,extra]=c.model.UpdObj(r.x,r.y,c.owner.ConPrb);
            end
            c.Ncal(level)=c.Ncal(level)+extra;c.g{level+1,1}=c.model.G;
            c.sampler.BayStc=c.model;
        end
        function c=select_run_seeds(~,c,r,level)
            if c.targets(level+1)<=0,error('SeMCS:ScheduleExhausted','Sample budget exhausted before seed selection.');end
            target=c.targets(level+1);[c.xs,c.ys]=SedSlt(c.model,r.x,c.y{level},target);
            c.nc=size(c.xs,2);c.steps=round(target/c.nc);
            if c.owner.Verbose
                fprintf('Level %d: N=%d, chains=%d, steps=%d, calls=%d\n',level,r.count,r.shape(1),r.shape(2),c.Ncal(level));
            end
        end
        function [yes,c]=run_complete(~,c,~,phase)
            yes=false;if ~strcmp(phase,'after_seeds'),return;end
            if c.nc<=0 || numel(c.x)==c.limit,yes=true;return;end
            c.complete=c.model.FlgCvg;yes=c.complete;
        end
        function out=finish_run(~,c)
            out=SeMCS.pack_output(c);
        end
    end
end
