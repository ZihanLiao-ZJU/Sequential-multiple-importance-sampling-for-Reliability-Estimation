classdef LSF_nD_LinearLSF2
    % Linear LSF with lognormal inputs (nD)
    % Configuration
    % ---------------------------------------------------------------------
    % 1. Dimension: n (user-specified)
    % 2. Input RVs: identically distributed lognormal (mu = 1.0, sigma = 0.2)
    % 3. LSF (in theta-space):
    %    G(theta) = (n + 3*sigma*sqrt(n)) - sum(theta)
    % 4. Failure event: G(theta) <= 0
    % 5. Pf_ref (varying with n, as reported):
    %    n = 10    : Pf_ref = 2.73e-3
    %    n = 20    : Pf_ref = 2.28e-3
    %    n = 50    : Pf_ref = 1.91e-3
    %    n = 100   : Pf_ref = 1.74e-3
    %    n = 1000  : Pf_ref = 1.47e-3
    % ---------------------------------------------------------------------
    % Reference:
    % ---------------------------------------------------------------------
    % [1] Dang, C., & Xu, J. Unified reliability assessment for problems with low- to
    %     high-dimensional random inputs using the Laplace transform and a mixture
    %     distribution. Reliability Engineering & System Safety, 204 (2020), 107124.
    %     https://doi.org/10.1016/j.ress.2020.107124
    % ---------------------------------------------------------------------

    properties(Constant)
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Nfun       length of function list                          double            [1,1]
        % -----------------------------------------------------------------------------------
        Nfun = 1;
    end

    properties
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Dist       Joint distribution (lognormal marginals + Nataf)  Nataf             [1,1]
        % Ndim       number of dimensions                             double            [1,1]
        % -----------------------------------------------------------------------------------
        Dist
        Ndim
    end

    methods
        function obj = LSF_nD_LinearLSF2(Ndim)
            obj.Ndim = Ndim;
            mu  = 1.0;
            sigma = 0.2;

            sLN2 = log(1 + (sigma^2)/(mu^2));
            sLN  = sqrt(sLN2);
            muLN = log(mu) - 0.5*sLN2;
            obj.Dist = repmat(makedist("Lognormal","mu",muLN,"sigma",sLN), obj.Ndim, 1);

            obj.Dist = Nataf(obj.Dist, eye(obj.Ndim));
        end

        function G = EvlLSF(obj, theta)
            % evaluate the function value
            % ----------------------------------------------------------------------
            % SYNTAX:
            % G = EvlLSF(obj,theta)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : constructed class
            % theta : variable samples                                   [Ndim,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % G     : function value; failure corresponds to G < 0                [1,Nsam]
            % ----------------------------------------------------------------------
            % REFERENCES:
            % ----------------------------------------------------------------------
            % [1] Dang & Xu (2020).
            % ----------------------------------------------------------------------
            sigma = 0.2;
            G = (obj.Ndim + 3*sigma*sqrt(obj.Ndim)) - sum(theta, 1);
        end
    end
end