classdef Ex6
    % Example 6: Brotone Bridge tower reliability problem (5D)   [paper Table 1]
    %
    % Definition
    % ---------------------------------------------------------------------
    % g(u) = X5 - 5*X3*X4*(47.4)^3 / ( 3*( 5*X1*X2 - 2*X3*(47.4)^2 ) )
    %           - 5*(1170*47.4)*X4*(0.4)^2 / ( 3*( 5*X1*X2 - 2*X3*(47.4)^2 ) )
    %           - 47.4*X4
    % Failure event:  g(u) < 0
    %
    % Configuration
    % ---------------------------------------------------------------------
    % 1) Dimension: 5
    % 2) Input RVs: independent with mixed normal + lognormal distributions (mean/COV specification)
    % 3) Reference failure probability: Pf_ref = 2.370e-4
    %
    % Reference
    % ---------------------------------------------------------------------
    % [43] C. Lan, H. Li, J. Peng, D. Sun, A structural reliability-based sensitivity analysis method using particle swarm optimization: relative convergence rate, Journal of Zhejiang University-SCIENCE A 17 (2016) 961-973. https://doi.org/10.1631/jzus.A1500255.
    %
    % Syntax
    % ---------------------------------------------------------------------
    % obj = Ex6
    properties
        Nfun = 1;
        Ndim = 5;
        Dist
    end

    methods
        function obj = Ex6()

            % -------------------------
            % Table 3: mean and COV
            % -------------------------
            mu  = [17.6; 4e7; 9.01e4; 2.12e3; 2.6e5];
            cov = [0.10; 0.08; 0.10;   0.15;   0.15];

            sigma = cov .* mu;

            dist = [];

            % X1 : Normal
            dist = [dist;
                makedist("Normal","mu",mu(1),"sigma",sigma(1))];

            % X2 : Lognormal (convert from mean/COV to log-space parameters)
            s2 = sqrt(log(1 + cov(2)^2));
            m2 = log(mu(2)) - 0.5*s2^2;
            dist = [dist;
                makedist("Lognormal","mu",m2,"sigma",s2)];

            % X3 : Normal
            dist = [dist;
                makedist("Normal","mu",mu(3),"sigma",sigma(3))];

            % X4 : Normal
            dist = [dist;
                makedist("Normal","mu",mu(4),"sigma",sigma(4))];

            % X5 : Normal
            dist = [dist;
                makedist("Normal","mu",mu(5),"sigma",sigma(5))];

            % independent
            obj.Dist = Nataf(dist, eye(obj.Ndim));
        end


        function g = EvlLSF(~, X)
            % Evaluate the performance function at samples X
            % --------------------------------------------------------------
            % INPUT:
            % X : samples in original variable space                    [5,   Nsam]
            % OUTPUT:
            % g : PF values; failure corresponds to g <= 0             [1,   Nsam]
            % --------------------------------------------------------------
            % Notes:
            % --------------------------------------------------------------
            % - This implementation follows the deterministic model parameters
            %   (H,l,W,h) defined below and evaluates g(X) directly in X-space.

            % deterministic parameters
            H = 70.5;
            l = 47.4;
            W = 1170;
            h = 0.4 * H;

            X1 = X(1,:);
            X2 = X(2,:);
            X3 = X(3,:);
            X4 = X(4,:);
            X5 = X(5,:);

            Den = 3 * (5*X1.*X2 - 2*X3*l^2);

            term1 = (5 * X3 .* X4 * l^3) ./ Den;
            term2 = (5 * W .* X4 * l * h^2) ./ Den;

            g = X5 - term1 - term2 - X4*l;
        end
    end
end
