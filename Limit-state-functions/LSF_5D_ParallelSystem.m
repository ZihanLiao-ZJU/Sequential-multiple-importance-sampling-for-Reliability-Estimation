classdef LSF_5D_ParallelSystem
    % 5-RV parallel system LSF in standard normal space
    % Configuration
    % ---------------------------------------------------------------------
    % 1. Dimension: 5
    % 2. Input RVs: independent standard normal
    % 3. LSF:
    %    G(u) = max{ 2.677 - u1 - u2,
    %               2.500 - u2 - u3,
    %               2.323 - u3 - u4,
    %               2.250 - u4 - u5 }
    % 4. Failure event: G(u) <= 0
    % 5. Pf_ref = 2.13e-4
    % ---------------------------------------------------------------------
    % Reference:
    % ---------------------------------------------------------------------
    % [30] Bjerager, P. Probability integration by directional simulation.
    %      Journal of Engineering Mechanics, 114 (1988), 1285–1302.
    %      https://doi.org/10.1061/(ASCE)0733-9399(1988)114:8(1285)
    % [31] Bourinet, J.-M. Rare-event probability estimation with adaptive support
    %      vector regression surrogates. Reliability Engineering & System Safety,
    %      150 (2016), 210–221. https://doi.org/10.1016/j.ress.2016.01.023.
    % ---------------------------------------------------------------------

    properties
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Nfun       length of function list                          double            [1,1]
        % Ndim       number of dimensions                             double            [1,1]
        % Dist       Joint distribution (independent standard normal)  Nataf             [1,1]
        % -----------------------------------------------------------------------------------
        Nfun = 1;
        Ndim = 5;
        Dist
    end

    methods
        function obj = LSF_5D_ParallelSystem()
            % Constructor: build independent standard normal marginals
            dist = repmat(makedist("Normal","mu",0,"sigma",1), obj.Ndim, 1);
            obj.Dist = Nataf(dist, eye(obj.Ndim));
        end

        function G = EvlLSF(~, u)
            % Evaluate the limit state function at samples u
            % ----------------------------------------------------------------------
            % SYNTAX:
            % G = EvlLSF(obj,u)
            % ----------------------------------------------------------------------
            % INPUTS:
            % u     : standard normal samples                           [Ndim,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % G     : LSF values; system failure corresponds to G <= 0   [1,Nsam]
            % ----------------------------------------------------------------------
            % Notes:
            % ----------------------------------------------------------------------
            % - Component failure is defined as gi <= 0 for each linear limit state.
            % - Parallel system: system fails only if ALL components fail, i.e.,
            %   max(gi) <= 0.

            % ---- component limit states ----
            g1 = 2.677 - u(1,:) - u(2,:);
            g2 = 2.500 - u(2,:) - u(3,:);
            g3 = 2.323 - u(3,:) - u(4,:);
            g4 = 2.250 - u(4,:) - u(5,:);

            % ---- parallel system aggregation ----
            G = max([g1; g2; g3; g4], [], 1);
        end
    end
end