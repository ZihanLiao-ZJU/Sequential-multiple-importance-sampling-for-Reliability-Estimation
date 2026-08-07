classdef Ex12
    % Example 12: linear PF with lognormal inputs (nD)   [paper Table 4]
    %
    % Definition
    % ---------------------------------------------------------------------
    % g(u) = ( n + 3*sigma*sqrt(n) ) - sum(u),   with lognormal inputs
    % X ~ Lognormal(mu = 1.0, sigma = 0.2) in physical space
    % Failure event:  g(u) < 0
    %
    % Configuration
    % ---------------------------------------------------------------------
    % 1) Dimension: n (user-specified)
    % 2) Input RVs: identically distributed lognormal (mu = 1.0, sigma = 0.2)
    % 3) Reference failure probability: Pf_ref = 1.91e-3 (n = 50; values for n = 10,...,1000 in the code)
    %
    % Reference
    % ---------------------------------------------------------------------
    % [49] C. Dang, J. Xu, Unified reliability assessment for problems with low- to high-dimensional random inputs using the Laplace transform and a mixture distribution, Reliability Engineering & System Safety 204 (2020) 107124. https://doi.org/10.1016/j.ress.2020.107124.
    %
    % Syntax
    % ---------------------------------------------------------------------
    % obj = Ex12(Ndim)
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
        function obj = Ex12(Ndim)
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
