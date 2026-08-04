classdef SeMIS_low_dimension
    % SeMIS with Nataf transformation for structural reliability problems
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
        function obj = SeMIS_low_dimension(LSF)
            % Initialize the SeMIS_low_dimension object
            obj.TF   = LSF;
            obj.Nfun = LSF.Nfun + 3;
            obj.Ndim = LSF.Ndim;

            beta   = 2.33;
            lambda = normpdf(beta) ./ (1 - normcdf(beta));
            sig    = sqrt((obj.Ndim + beta * lambda) ./ obj.Ndim);

            obj.G.p      = 0.1;
            obj.G.Ite    = 0;
            obj.G.fth    = inf;
            obj.G.Krg    = {};
            obj.G.krg_th = [];
            obj.G.kappa  = 2.0;

            obj.G.L  = eye(obj.Ndim) * sig;
            obj.G.CV = obj.G.L * obj.G.L.';
            obj.G.mu = zeros(obj.Ndim,1);

            obj.Ncal = 0;
        end

        function obj = EvlY(obj,u)
            % Evaluate the function values associated with the input samples
            u = u(1:end,:);
            obj.X = u;
            N_sam = size(u,2);

            % Cheap component used in the standard normal space
            logp = logGauss(u);

            % Screening condition defined by kriging in the U-space
            if obj.G.Ite > 0
                ind = true(1,N_sam);
                for i = 1:obj.G.Ite
                    [mu_pred,sd_pred] = predict(obj.G.Krg{i},u.');
                    score = (mu_pred - obj.G.kappa .* sd_pred).';
                    ind = ind & (score <= obj.G.krg_th(i));
                end
            else
                ind = true(1,N_sam);
            end

            % Default setting: no LSF evaluation is performed
            f    = inf(1,N_sam);
            ncal = 0;

            if any(ind)
                x      = U2X(obj,u(:,ind));
                f(ind) = obj.TF.EvlLSF(x);
                ncal   = sum(ind);
            end

            obj.Y    = [logp; f];
            obj.Ncal = ncal;
        end

        function Li = EvlLKF(obj)
            % Evaluate the intermediate log-likelihood
            f     = obj.Y(2,:);
            fth   = obj.G.fth;
            N_sam = size(f,2);

            if isinf(fth)
                Li = zeros(1,N_sam);
                return;
            end

            if obj.G.Ite > 0
                ind = true(1,N_sam);
                for i = 1:obj.G.Ite
                    [mu_pred,sd_pred] = predict(obj.G.Krg{i},obj.X.');
                    score = (mu_pred - obj.G.kappa .* sd_pred).';
                    ind = ind & (score <= obj.G.krg_th(i));
                end
            else
                ind = true(1,N_sam);
            end

            Li = log(double((f < fth) & ind));
        end

        function Pi = EvlPDF(obj)
            % Evaluate the intermediate log-prior density pi
            Pi = logGauss(obj.X,obj.G.mu,obj.G.CV);
        end

        function [y,g,obj,Ncal] = UpdObj(obj,u,y)
            % Update the intermediate distribution parameters
            f = y(4,:);
            u = u(1:end,:);
            g = obj.G;

            fth = max(quantile(f,g.p),0);
            ind_sed = f < fth;

            % Fit kriging using all truly evaluated samples
            ind_fit = isfinite(f);
            u_fit   = u(:,ind_fit);
            f_fit   = f(ind_fit);

            mdl = fitrgp( ...
                u_fit.',f_fit.', ...
                'KernelFunction','ardsquaredexponential', ...
                'BasisFunction','constant', ...
                'Standardize',true, ...
                'FitMethod','exact', ...
                'PredictMethod','exact');

            % Conservative kriging score
            [mu_pred,sd_pred] = predict(mdl,u.');
            score = (mu_pred - g.kappa .* sd_pred).';

            % Threshold chosen so that all current seeds are retained
            krg_th_tmp = max(score(ind_sed));

            g.Krg{end+1} = mdl;
            g.krg_th     = cat(1,g.krg_th,krg_th_tmp);
            g.fth        = fth;
            g.Ite = g.Ite + 1;

            obj.G   = g;
            Ncal    = 0;
        end

        function FlgCvg = get.FlgCvg(obj)
            % Check whether the convergence criterion is satisfied
            FlgCvg = obj.G.fth <= 0;
        end
    end
end