classdef LSF_3D_Benchmarkp6
    % Polynomial benchmark LSF (3D) in standard normal space
    % Configuration
    % ---------------------------------------------------------------------
    % 1. Dimension: 3
    % 2. Pf_ref = 3.13e-4
    % ---------------------------------------------------------------------
    % Reference:
    % ---------------------------------------------------------------------
    % [1] Guan, X. L., & Melchers, R. E. Effect of response surface parameter variation on
    %     structural reliability estimates. Structural Safety, 23 (2001), 429–444.
    %     https://doi.org/10.1016/S0167-4730(02)00013-9
    % ---------------------------------------------------------------------
    properties
        % parameters and model settings
        % ---------------------------------------------------------------------
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Nfun       length of function list                          double            [1,1]
        % Ndim       number of dimensions                             double            [1,1]
        % Dist       Reference distribution (standard normal space)    Nataf             [1,1]
        % -----------------------------------------------------------------------------------
        Nfun = 1;
        Ndim = 3;
        Dist = Nataf(repmat(makedist("Normal","mu",0,"sigma",1),3,1),eye(3));
    end

    methods
        function L = EvlLSF(~,theta)
            % Evaluate the limit state function at samples theta
            % ----------------------------------------------------------------------
            % SYNTAX:
            % L = obj.EvlLSF(theta)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : constructed class
            % theta : samples in standard normal space                        [Ndim,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % L     : LSF values; failure corresponds to L < 0                [1,Nsam]
            % ----------------------------------------------------------------------
            % REFERENCES:
            % ----------------------------------------------------------------------
            % [1] Guan, X. L., & Melchers, R. E. (2001).
            % ----------------------------------------------------------------------

            % evaluation
            L = 1/40 *theta(1,:).^4 + 2*theta(2,:).^2 + theta(3,:)+3;
        end

       
    end
end