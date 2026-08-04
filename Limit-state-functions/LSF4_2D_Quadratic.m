classdef LSF4_2D_Quadratic
    % A 2D example with four LSFs
    % Configuration
    % ---------------------------------------------------------------------
    % 1. Dimension: 2
    % 2. Pf_1 = 8.53e-4;
    %    Pf_2 = 2.29e-4;
    %    Pf_3 = 8.79e-4;
    %    Pf_4 = 2.40e-4;
    % ---------------------------------------------------------------------
    % Reference:
    % ---------------------------------------------------------------------
    % [1] Li, Hong-Shuang, Yuan-Zhuo Ma, and Zijun Cao. "A generalized subset simulation approach for estimating small failure probabilities of multiple stochastic responses." Computers & Structures 153 (2015): 239-251.
    %     -- Section 4.2
    % ---------------------------------------------------------------------
    properties
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Parm       parameters                                       double         [Npar,1]
        % Nfun       length of function list                          double            [1,1]
        % Ndim       number of dimensions                             double      [Ncha,Ndat]
        % Rdis       Reference distribution                           /                 [1,1]
        % -----------------------------------------------------------------------------------
        Nfun = 4;
        Ndim = 2;
        Dist = Nataf([makedist("Normal","mu",0,"sigma",1);makedist("Normal","mu",0,"sigma",1)],eye(2));
    end

    methods
        function L = EvlLSF(~,theta)
            % Evaluate the function value
            % ----------------------------------------------------------------------
            % SYNTAX:
            % L = EvlFun(obj,theta)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : constructed class
            % theta : variable samples                                   [Ndim,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % L     : LSF value                                          [Nfun,Nsam]
            % ----------------------------------------------------------------------

            % initialization
            L1 = 3+(theta(1,:)-theta(2,:)).^2/10-(theta(1,:)+theta(2,:))/sqrt(2);
            L2 = theta(1,:)-theta(2,:)+7/sqrt(2);
            L3 = 3+(theta(1,:)-theta(2,:)).^2/10+(theta(1,:)+theta(2,:))/sqrt(2);
            L4 = theta(2,:)-theta(1,:)+7/sqrt(2);
            L = [L1;L2;L3;L4];
        end
    end
end
