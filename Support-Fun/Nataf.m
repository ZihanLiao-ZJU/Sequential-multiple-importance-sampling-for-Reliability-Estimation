classdef Nataf
    % Nataf transformation
    properties
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Marg_X     Marginal distribution of input variables X       class          [Ndim,1]
        % Rho_X      Correlation matrix of standard normal            double      [Ndim,Ndim]
        %            variables X
        % Rho_U      Correlation matrix of standard normal            double      [Ndim,Ndim]
        %            variables U
        % RhoL_U     Lower triangular matrix of RhoL_U                double      [Ndim,Ndim]
        % RhoL_X     Lower triangular matrix of RhoL_X                double      [Ndim,Ndim]
        % -----------------------------------------------------------------------------------
        Marg_X
        Rho_X
        Rho_U
        RhoL_U
        RhoL_X
        Ndim
    end

    properties(Constant)
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Nfun       Length of function list                          double           [1,1]
        % LogP       output of PDF with logrithm (default: true)      logical           [1,1]
        % -----------------------------------------------------------------------------------
        Nfun = 1;
    end

    methods

        function obj = Nataf(Marg_X,Rho_X)
            % Initilization
            % -------------------------------------------------------------------------------
            % SYNTAX:
            % obj = Nataf(marg_x,Rho_X)
            % -------------------------------------------------------------------------------
            % INPUTS:
            % marg_x  : Marginal distribution of input variables x
            % -------------------------------------------------------------------------------
            % OUTPUTS:
            % obj     : constructed class
            % -------------------------------------------------------------------------------

            % extraction from inputs
            obj.Marg_X = Marg_X;
            obj.Rho_X = Rho_X;
            Ndim = size(Marg_X,1);
            obj.Ndim = Ndim;
            % evaluate RhoL_X
            % --------------------------------------------------
            % perform Choleski decomposition
            [obj.RhoL_X,flg_pos] = chol(Rho_X,'lower');
            if flg_pos ~= 0
                error('Correlation matrix of variables is not positive definite');
            end
            % --------------------------------------------------
            % evaluate RhoL_U
            % --------------------------------------------------
            % initilization of RhoL_U
            if all(Rho_X==eye(Ndim))
                Rho_U = eye(Ndim);
            else
                Rho_U = eye(Ndim);
                for i = 1:Ndim
                    for j = i+1:Ndim
                        % explicit expression of correlatipn transformation between lognormal and normal
                        if strcmpi(Marg_X(i).DistributionName,'Normal') && strcmpi(Marg_X(j).DistributionName,'Normal')
                            Rho_U(i,j) = Rho_X(i,j);
                            Rho_U(j,i) = Rho_X(i,j);
                        elseif strcmpi(Marg_X(i).DistributionName,'Normal') && strcmpi(Marg_X(j).DistributionName,'Lognormal')
                            cv_j = Marg_X(j).std/Marg_X(j).mean;
                            Rho_U(i,j) = Rho_X(i,j)*cv_j/sqrt(log(1+cv_j^2));
                            Rho_U(j,i) = Rho_U(i,j);
                        elseif strcmpi(Marg_X(i).DistributionName,'Lognormal') && strcmpi(Marg_X(j).DistributionName,'Normal')
                            cv_i = Marg_X(i).std/Marg_X(i).mean;
                            Rho_U(i,j) = Rho_X(i,j)*cv_i/sqrt(log(1+cv_i^2));
                            Rho_U(j,i) = Rho_U(i,j);
                        elseif strcmpi(Marg_X(i).DistributionName,'Lognormal') && strcmpi(Marg_X(j).DistributionName,'Lognormal')
                            cv_i = Marg_X(i).std/Marg_X(i).mean;
                            cv_j = Marg_X(j).std/Marg_X(j).mean;
                            Rho_U(i,j) = log(1+Rho_X(i,j)*cv_i*cv_j)/sqrt(log(1+cv_i^2)*log(1+cv_j^2));
                            Rho_U(j,i) = Rho_U(i,j);
                        else
                            % computes: the error of correlation transformation in Nataf
                            marg_x = [Marg_X(i);Marg_X(j)];
                            rho_x = Rho_X(i,j);

                            % Calculation of the transformed correlation matrix. This is achieved by a
                            % quadratic two-dimensional Gauss-Legendre integration
                            opt_GIP.u_max = [8;8];    % Integration bounds
                            opt_GIP.u_min = [-8;-8];   % Integration bounds
                            Nint    = [200;200];    % number of integration points along each dimension

                            [u,w] = GauIntPot(Nint,opt_GIP);
                            u = u';
                            w = w';

                            % Transformation of parameters of the distributions for the
                            % matlab cdf/icdf functions and calculations of those
                            fi = (marg_x(1).icdf(normcdf(u(:,1)))-marg_x(1).mean)/marg_x(1).std;
                            fj = (marg_x(2).icdf(normcdf(u(:,2)))-marg_x(2).mean)/marg_x(2).std;
                            coef = fi.*fj.*w;

                            rhoU2X_err  = @(rho_u)abs(rho_x-sum(coef./(2*pi*sqrt(1-rho_u.^2)).*exp( -1./(2*(1-rho_u.^2)).*(u(:,1).^2 - 2*rho_u.*u(:,1).*u(:,2) + u(:,2).^2) ),1));
                            Rho_U(i,j) = fminbnd(rhoU2X_err,0,1);
                            Rho_U(j,i) = Rho_U(i,j);
                        end
                    end
                end
            end
            Rho_U = max(Rho_U,-1);
            Rho_U = min(Rho_U,1);
            obj.Rho_U = Rho_U;
            % perform Choleski decomposition
            [obj.RhoL_U,flg_pos] = chol(Rho_U,'lower');
            if flg_pos ~= 0
                error('Correlation matrix of standard normal variables is not positive definite');
            end
            % --------------------------------------------------
        end

        function X = U2X(obj,U)
            % This function computes: standard normal space ===> variable sapce
            % Syntax
            % -------------------------------------------------------------------------------
            % X = U2X(obj,U)
            % -------------------------------------------------------------------------------
            %
            % Inputs:
            % -------------------------------------------------------------------------------
            % 1. obj     : constructed nataf class
            % 2. U       : standard normal samples                                   [Ndim,:]
            % -------------------------------------------------------------------------------
            %
            % Outputs:
            % -------------------------------------------------------------------------------
            % 1. X       : samples conform to the input distribution                 [Ndim,:]
            % -------------------------------------------------------------------------------

            SzeU = size(U);
            % Correlated standard normal variables
            Z = pagemtimes(obj.RhoL_U, U);
            X = zeros(SzeU);

            for idim = 1:obj.Ndim
                % Distribution name
                distName = obj.Marg_X(idim).DistributionName;
                if strcmpi(distName,'Normal') || strcmpi(distName,'正态')
                    % Normal distribution: linear transformation
                    mu    = obj.Marg_X(idim).mu;
                    sigma = obj.Marg_X(idim).sigma;
                    X(idim,:) = mu + sigma .* Z(idim,:);

                elseif strcmpi(distName,'Lognormal') || strcmpi(distName,'对数正态')
                    % Lognormal distribution: explicit variable transformation
                    % X = exp(mu_ln + sigma_ln * Z)
                    mu_ln    = obj.Marg_X(idim).mu;
                    sigma_ln = obj.Marg_X(idim).sigma;
                    X(idim,:) = exp(mu_ln + sigma_ln .* Z(idim,:));

                elseif strcmpi(distName,'Generalized Extreme Value') || strcmpi(distName,'广义极值')
                    % Generalized Extreme Value distribution: explicit inverse transform
                    % p = Phi(Z)
                    % If k == 0: Gumbel for maxima (right-skewed)
                    %   X = mu - sigma * log(-log(p))
                    % If k ~= 0:
                    %   X = mu + (sigma/k) * ( (-log(p)).^(-k) - 1 )
                    p = normcdf(Z(idim,:));
                    p = min(max(p, 1e-15), 1-1e-15);

                    k     = obj.Marg_X(idim).k;
                    sigma = obj.Marg_X(idim).sigma;
                    mu    = obj.Marg_X(idim).mu;

                    if abs(k) < 1e-12
                        % Gumbel (k = 0): maxima
                        X(idim,:) = mu - sigma .* log(-log(p));
                    else
                        % General GEV (k ~= 0)
                        X(idim,:) = mu + (sigma./k) .* ( (-log(p)).^(-k) - 1 );
                    end

                else
                    % General distribution: inverse CDF mapping
                    X(idim,:) = obj.Marg_X(idim).icdf( normcdf(Z(idim,:)) );
                end
            end
        end

        function U = X2U(obj,X)
            % This function computes: variable sapce ===> standard normal space
            % Syntax
            % -------------------------------------------------------------------------------
            % U = X2U(obj,X)
            % U = X2U(obj,X,opt_X2U)
            % -------------------------------------------------------------------------------
            %
            % Inputs:
            % -------------------------------------------------------------------------------
            % 1. obj     : constructed nataf class
            % 2. X       : samples conform to the input distribution              [Ndim,Nsam]
            % 3. opt_X2U : optional inputs
            % -------------------------------------------------------------------------------
            %
            % Outputs:
            % -------------------------------------------------------------------------------
            % 1. U       : standard normal samples                                [Ndim,Nsam]
            % -------------------------------------------------------------------------------

            SzeX = size(X);
            Z = zeros(SzeX);
            for idim = 1:obj.Ndim
                Z(idim,:) = norminv(obj.Marg_X(idim).cdf(X(idim,:)));
            end
            U = pagemtimes(inv(obj.RhoL_U),Z);
        end

        function pdv = EvlPDF(obj,X)
            % This function computes the (log) pdf of input samples
            % Syntax
            % -------------------------------------------------------------------------------
            % pdv = EvlPDF(obj,X)
            % -------------------------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------------------------
            % obj     : constructed nataf class
            % X       : samples conform to the input distribution                 [Ndim,Nsam]
            % -------------------------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------------------------
            % pdv     : joint (log) pdf value of X                                   [1,Nsam]
            % -------------------------------------------------------------------------------

            ndim = obj.Ndim;
            Nsam = size(X,2);
            % standard normal with the same marginal CDF to X
            Z = zeros(ndim,Nsam);
            % pdf of Z
            pdv_Z = zeros(ndim,Nsam);
            % pdf of Y(pdf of X with independent joint distribution)
            pdv_Y = zeros(ndim,Nsam);

            for idim = 1:ndim
                Z(idim,:) = norminv(obj.Marg_X(idim).cdf(X(idim,:)));
                pdv_Z(idim,:) = logGauss(Z(idim,:));
                pdv_Y(idim,:) = log(obj.Marg_X(idim).pdf(X(idim,:)));
            end
            % pdf pf U (pdf of Z multivariate with mean 0 std 1 and Rho_U)
            pdv_U = logGauss(Z,zeros(ndim,1),obj.Rho_U);
            pdv = sum(pdv_Y,1)-sum(pdv_Z,1)+pdv_U;
        end

        function cdv = EvlCDF(obj,X)
            % This function computes the (log) pdf of input samples
            % Syntax
            % -------------------------------------------------------------------------------
            % pdv = EvlCDF(obj,X)
            % -------------------------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------------------------
            % 1. obj     : constructed nataf class
            % 2. X       : samples conform to the input distribution              [Ndim,Nsam]
            % 3. logP    : pdf in logarithm or not (default:false)
            % -------------------------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------------------------
            % 1. cdv     : joint (log) cdf value of X                                [1,Nsam]
            % -------------------------------------------------------------------------------

            ndim = obj.Ndim;
            Nsam = size(X,2);
            % standard normal with the same marginal CDF to X
            Z = zeros(ndim,Nsam);
            for idim = 1:ndim
                Z(idim,:) = norminv(obj.Marg_X(idim).cdf(X(idim,:)));
            end
            cdv = logGauss(Z,zeros(ndim,1),obj.Rho_U);
        end

        function X = SamGen(obj,N)
            % Generate random samples with input size
            % Syntax
            % -------------------------------------------------------------------------------
            % X = SamGen(obj,N)
            % -------------------------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------------------------
            % obj  : constructed nataf class
            % N    : number of generated samples
            % -------------------------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------------------------
            % X    : generated samples                                            [Ndim,N]
            % -------------------------------------------------------------------------------
            % generation of randn numbers
            ndim = obj.Ndim;
            u = randn(ndim,N);
            X = obj.U2X(u);
        end
    end
end
