classdef SeMIS_high_dimension
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
        function obj = SeMIS_high_dimension(LSF)
            % Initialize the SeMIS_SuS object
            obj.TF  = LSF;
            obj.Nfun = LSF.Nfun + 3;
            obj.Ndim = LSF.Ndim;

            beta = 2.33;
            lambda = normpdf(beta) ./ (1-normcdf(beta));
            sig  = sqrt((obj.Ndim+(beta*lambda))./obj.Ndim);

            obj.G.p    = 0.1;
            obj.G.Ite    = 0;
            obj.G.fth  = inf;
            obj.G.Beta = [];
            obj.G.f_hat_th = [];

            obj.G.L  = eye(obj.Ndim) * sig;
            obj.G.CV = obj.G.L*obj.G.L.';
            obj.G.mu = zeros(obj.Ndim,1);

            obj.Ncal = 0;
        end

        function obj = EvlY(obj,u)
            % Evaluate the function values associated with the input samples
            u = u(1:end,:);
            obj.X = u;
            N_sam = size(u,2);

            % Cheap component used to match the Gaussian proposal N(mu,CV)
            logp = logGauss(u);

            % Screening condition defined in the U-space
            if obj.G.Ite > 0
                f_hat   = obj.G.Beta.' * [ones(1,N_sam); u];       % [1,Nsam]
                ind     = all(f_hat <= obj.G.f_hat_th,1);
            else
                ind = true(1,N_sam);
            end

            % Default setting: no LSF evaluation is performed
            f    = inf(1,N_sam);
            ncal = 0;

            if any(ind)
                x      = U2X(obj, u(:,ind));
                f(ind) = obj.TF.EvlLSF(x);
                ncal   = sum(ind);
            end

            obj.Y    = [logp; f];
            obj.Ncal = ncal;
        end

        function Li = EvlLKF(obj)
            % Evaluate the intermediate log-likelihood
            f    = obj.Y(2,:);
            fth = obj.G.fth;
            N_sam = size(f,2);

            if isinf(fth)
                Li = zeros(1, N_sam);
                return;
            end

            f_hat = obj.G.Beta.' * [ones(1,N_sam); obj.X];             % [1,Nsam]
            Li    = log(double((f < fth) & (all(f_hat <= obj.G.f_hat_th,1))));
        end

        function Pi = EvlPDF(obj)
            % Evaluate the intermediate log-prior density pi
            Pi = logGauss(obj.X, obj.G.mu, obj.G.CV);
        end

        function [y,g,obj,Ncal] = UpdObj(obj,u,y)
            % Update the intermediate distribution parameters
            f = y(4,:);
            u = u(1:end,:);
            g = obj.G;
            N_sam = size(u,2);

            fth = max(quantile(f, g.p), 0);
            ind_sed = f < fth;
            
            u_sed = u(:,ind_sed);
            f_sed = f(:,ind_sed);
            [~,~,~,~,Beta_tmp] = plsregress(u_sed.', f_sed.', 1);
            f_hat  = Beta_tmp.' * [ones(1,N_sam); u];
            f_hat_tmp =  max(f_hat(:,ind_sed),[],2);

            g.Beta = cat(2,g.Beta,Beta_tmp);
            g.f_hat_th  = cat(1,g.f_hat_th,f_hat_tmp);
            g.fth  = fth;
            g.Ite = g.Ite + 1;

            obj.G   = g;
            Ncal    = 0;
        end

        function FlgCvg = get.FlgCvg(obj)
            % Check whether the convergence criterion is satisfied
            FlgCvg = obj.G.fth<=0;
        end
    end
end