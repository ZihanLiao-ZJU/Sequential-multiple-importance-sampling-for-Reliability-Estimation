classdef LSF_7D_Reinforced_concrete_beam
    % Reinforced concrete beam reliability problem (7D)
    % Configuration
    % ---------------------------------------------------------------------
    % 1. Dimension: 7
    % 2. Independent random variables
    % 3. Mixed Gumbel (GEV with k=0) + normal + lognormal distributions
    % 4. LSF:
    %    Z = X2*X3*X4 - (X5*X3^2*X4^2)/(X6*X7) - X1
    % 5. Failure event: Z <= 0
    % 6. Pf_ref = 3.20e-3
    % ---------------------------------------------------------------------
    % Reference:
    % ---------------------------------------------------------------------
    % [1] Gong, J., Yi, P., & Zhao, N. Non-gradient–based algorithm for structural
    %     reliability analysis. Journal of Engineering Mechanics, 140 (2014), 04014029.
    %     https://doi.org/10.1061/(ASCE)EM.1943-7889.0000722
    % ---------------------------------------------------------------------

    properties
        Nfun = 1;
        Ndim = 7;
        Dist
    end

    methods
        function obj = LSF_7D_Reinforced_concrete_beam()

            dist = [];

            % -------------------------
            % X1 : Gumbel
            % -------------------------
            mu_mean = 0.01;
            sig_std = 0.003;

            gamma = 0.5772156649015329;
            beta  = sig_std * sqrt(6) / pi;
            mu_gev = mu_mean - gamma * beta;

            dist = [dist;makedist('GeneralizedExtremeValue','k',0,'sigma',beta,'mu',mu_gev)];

            % -------------------------
            % X2 : Normal
            % -------------------------
            dist = [dist;
                makedist("Normal","mu",0.30,"sigma",0.015)];

            % -------------------------
            % X3 : Normal
            % -------------------------
            dist = [dist;
                makedist("Normal","mu",360.0,"sigma",36.0)];

            % -------------------------
            % X4 : Lognormal
            % -------------------------
            mu4 = 226e-6;
            sd4 = 11.3e-6;

            s4 = sqrt(log(1 + (sd4/mu4)^2));
            m4 = log(mu4) - 0.5*s4^2;

            dist = [dist;
                makedist("Lognormal","mu",m4,"sigma",s4)];

            % -------------------------
            % X5 : Normal
            % -------------------------
            dist = [dist;
                makedist("Normal","mu",0.5,"sigma",0.05)];

            % -------------------------
            % X6 : Normal
            % -------------------------
            dist = [dist;
                makedist("Normal","mu",0.12,"sigma",0.006)];

            % -------------------------
            % X7 : Lognormal
            % -------------------------
            mu7 = 40;
            sd7 = 6;

            s7 = sqrt(log(1 + (sd7/mu7)^2));
            m7 = log(mu7) - 0.5*s7^2;

            dist = [dist;
                makedist("Lognormal","mu",m7,"sigma",s7)];

            obj.Dist = Nataf(dist, eye(obj.Ndim));   % independent
        end


        function Z = EvlLSF(~, X)
            % Evaluate the limit state function at samples X
            % ---------------------------------------------------------------------
            % INPUT:
            % X : samples in original variable space                         [7,   Nsam]
            % OUTPUT:
            % Z : LSF values; failure corresponds to Z <= 0                  [1,   Nsam]
            % ---------------------------------------------------------------------

            X1 = X(1,:);
            X2 = X(2,:);
            X3 = X(3,:);
            X4 = X(4,:);
            X5 = X(5,:);
            X6 = X(6,:);
            X7 = X(7,:);

            Den = X6 .* X7;

            Z = X2.*X3.*X4 - (X5 .* X3.^2 .* X4.^2) ./ Den - X1;
        end
    end
end