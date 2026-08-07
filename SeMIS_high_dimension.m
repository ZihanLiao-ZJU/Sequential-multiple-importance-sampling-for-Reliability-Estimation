classdef SeMIS_high_dimension
    % SeMIS with Nataf transformation for structural reliability problems
    % High-dimensional variant: PLS-based hybrid truncation (Eqn. 26).
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
        function obj = SeMIS_high_dimension(PF)
            % Initialize the SeMIS_high_dimension object
            obj.TF   = PF;
            obj.Nfun = PF.Nfun + 3;
            obj.Ndim = PF.Ndim;

            % Approximate optimal variance-inflation factor (Eqn. 19):
            %   sigma_tilde* = sqrt(1 + beta * lambda(beta) / n)
            beta      = 2.33;
            lambda_beta = normpdf(beta) ./ (1 - normcdf(beta));
            sigma     = sqrt((obj.Ndim + beta * lambda_beta) ./ obj.Ndim);

            obj.G.p      = 0.1;      % intermediate level probability p (Section 4.1.2)
            obj.G.Ite    = 0;        % current iteration counter
            obj.G.l      = inf;      % exact PF threshold l_i (Eqn. 20)
            obj.G.pls_beta = [];     % PLS coefficients [b; a] (Eqn. 26), one column per iteration
            obj.G.l_hat    = [];     % surrogate thresholds l_hat_m (Eqn. 27)

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

            % Screening condition defined in the U-space by the PLS hyperplanes:
            %   g_hat_m(u) = a_m' * u + b_m  (Eqn. 26)
            if obj.G.Ite > 0
                g_hat   = obj.G.pls_beta.' * [ones(1,n_sam); u];  % [1,Nsam]
                ind_scr = all(g_hat <= obj.G.l_hat,1);
            else
                ind_scr = true(1,n_sam);
            end

            % Default setting: no PF evaluation is performed
            g_val = inf(1,n_sam);
            n_cal = 0;

            if any(ind_scr)
                x      = U2X(obj, u(:,ind_scr));
                g_val(ind_scr) = obj.TF.EvlLSF(x);
                n_cal  = sum(ind_scr);
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
                log_L = zeros(1, n_sam);
                return;
            end

            g_hat = obj.G.pls_beta.' * [ones(1,n_sam); obj.X];   % PLS prediction (Eqn. 26)
            log_L = log(double((g_val < l) & (all(g_hat <= obj.G.l_hat,1))));
        end

        function log_pi = EvlPDF(obj)
            % Evaluate the intermediate log-prior density = log rPDF pi_i(u) (Eqn. 12)
            log_pi = logGauss(obj.X, obj.G.mu_pi, obj.G.Sigma);
        end

        function [y,g,obj,N_cal] = UpdObj(obj,u,y)
            % Update the intermediate distribution parameters
            %   l          : exact PF threshold = p-quantile of g (Eqn. 24)
            %   pls_beta   : PLS hyperplane fitted on the seed samples only (Eqn. 26)
            %   l_hat      : conservative surrogate threshold = max over seeds (Eqn. 27)
            g_val = y(4,:);
            u = u(1:end,:);
            g = obj.G;
            n_sam = size(u,2);

            l        = max(quantile(g_val, g.p), 0);   % threshold l_i (Eqn. 24)
            ind_seed = g_val < l;                      % seed set S_{i+1} (Eqn. 24)

            u_seed = u(:,ind_seed);
            g_seed = g_val(:,ind_seed);
            [~,~,~,~,pls_beta] = plsregress(u_seed.', g_seed.', 1);   % PLS coefficients (Eqn. 26)
            g_hat  = pls_beta.' * [ones(1,n_sam); u];
            l_hat_tmp = max(g_hat(:,ind_seed),[],2);   % conservative threshold (Eqn. 27)

            g.pls_beta = cat(2,g.pls_beta,pls_beta);
            g.l_hat    = cat(1,g.l_hat,l_hat_tmp);
            g.l        = l;
            g.Ite = g.Ite + 1;

            obj.G   = g;
            N_cal   = 0;
        end

        function FlgCvg = get.FlgCvg(obj)
            % Check whether the convergence criterion is satisfied (Section 4.1.2)
            FlgCvg = obj.G.l <= 0;
        end
    end
end
