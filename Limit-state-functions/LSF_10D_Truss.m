classdef LSF_10D_Truss
    % 10-bar truss reliability problem (10D)
    % Configuration
    % ---------------------------------------------------------------------
    % 1. Dimension: 10
    % 2. Input RVs (independent via Nataf):
    %    - Loads F1..F6: Gumbel (GEV with k = 0), mean = 50 kN, std = 7.5 kN
    %    - Young's modulus E1..E2: lognormal, mean = 210 GPa, std = 21 GPa
    %    - Areas A1..A2: lognormal, (20e-4, 2e-4) and (10e-4, 1e-4)
    % 3. LSF:
    %    L = 0.14 - |u(20)|
    %    where u(20) is the target displacement component from the FEM solver
    % 4. Failure event: L <= 0
    % 5. Pf_ref = 3.45e-5
    % ---------------------------------------------------------------------
    % Reference:
    % ---------------------------------------------------------------------
    % [39] Sun, Z., Wang, J., Li, R., & Tong, C. LIF: A new Kriging based learning
    %      function and its application to structural reliability analysis.
    %      Reliability Engineering & System Safety, 157 (2017), 152–165.
    %      https://doi.org/10.1016/j.ress.2016.09.003
    % [40] Schuëller, G. I., & Pradlwarter, H. J. Benchmark study on reliability
    %      estimation in higher dimensions of structural systems – An overview.
    %      Structural Safety, 29 (2007), 167–182.
    %      https://doi.org/10.1016/j.strusafe.2006.07.010
    % ---------------------------------------------------------------------

    properties(Constant)
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Nfun       length of function list                          double            [1,1]
        % Ndim       number of dimensions                             double            [1,1]
        % -----------------------------------------------------------------------------------
        Nfun = 1;
        Ndim = 10;
    end

    properties
        % properties:
        % Name       Description                                      Type               Size
        % -----------------------------------------------------------------------------------
        % Dist       Joint distribution (mixed marginals + Nataf)      Nataf             [1,1]
        % Fem        FEM model template (deterministic topology)       FEM              [1,1]
        % -----------------------------------------------------------------------------------
        Dist
        Fem
    end

    methods
        function obj = LSF_10D_Truss
            % distribution of load
            mu = 50e3;
            sig = 7.5e3;
            gamma_euler = 0.57721566490153286060;
            % Calculate beta using the standard deviation
            beta = sig * sqrt(6) / pi;
            % Calculate mode (mu) using the mean and beta
            mode = mu - beta * gamma_euler;
            dist_F = repmat(makedist('GeneralizedExtremeValue','k',0,'sigma',beta,'mu',mode),6,1);

            % distribution of E
            mu = 210e9;
            sig = 21e9;
            cov = sig/mu;
            sig_log = log(1+cov^2)^0.5;
            mode = log(mu) - 0.5*sig_log^2;
            dist_E = repmat(makedist('Lognormal','mu',mode,'sigma',sig_log),2,1);

            % distribution of A1
            mu = 20e-4;
            sig = 2e-4;
            cov = sig/mu;
            sig_log = log(1+cov^2)^0.5;
            mode = log(mu) - 0.5*sig_log^2;
            dist_A1 = makedist('Lognormal','mu',mode,'sigma',sig_log);

            % distribution of A2
            mu = 10e-4;
            sig = 1e-4;
            cov = sig/mu;
            sig_log = log(1+cov^2)^0.5;
            mode = log(mu) - 0.5*sig_log^2;
            dist_A2 = makedist('Lognormal','mu',mode,'sigma',sig_log);

            obj.Dist = [dist_F;dist_E;dist_A1;dist_A2];
            obj.Dist = Nataf(obj.Dist,eye(10));

            % FEM construction
            % Nodes
            Node = zeros(13,3);
            for i = 1:6
                Node(i,:) = [2+4*(i-1),2,0];
            end
            for i = 7:13
                Node(i,:) = [4*(i-7),0,0];
            end
            % Elements
            Ele = zeros(23,2);
            Ele_pro = repmat([1,1,20e-4,210e9],23,1);
            for i=1:5
                Ele(i,:) = [i,i+1];
            end
            for i=6:11
                Ele(i,:) = [i+1,i+2];
            end
            for i = 12:17
                Ele(i,:) = [i-11,i-5];
            end
            for i = 18:23
                Ele(i,:) = [i-17,i-10];
            end
            % Boundary condiitons
            BC = [7,1,0
                7,2,0
                13,2,0];

            % load
            F = zeros(6,3);
            for i = 1:6
                F(i,:) = [i,2,-50e3];
            end

            obj.Fem = FEM(Node,Ele,Ele_pro,BC,F);
        end

        function L = EvlLSF(obj,theta)
            % evaluate the function value
            % ----------------------------------------------------------------------
            % SYNTAX:
            % L = EvlLSF(obj,theta)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : constructed class
            % theta : variable samples                                   [Ndim,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % L     : function value; failure corresponds to L < 0                [1,Nsam]
            % ----------------------------------------------------------------------
            % REFERENCES:
            % ----------------------------------------------------------------------
            % [39] Sun et al. (2017); [40] Schuëller & Pradlwarter (2007).
            % ----------------------------------------------------------------------

            Nsam = size(theta,2);
            Ele_pro = obj.Fem.Ele_pro;
            Node = obj.Fem.Node;
            Ele = obj.Fem.Ele;
            BC = obj.Fem.BC;

            L = zeros(1,Nsam);
            for isam = 1:Nsam
                Ele_pro_tmp = Ele_pro;
                theta_tmp = theta(:,isam);
                F_sam = theta_tmp(1:6);
                E_sam = theta_tmp(7:8);
                A_sam = theta_tmp(9:10);

                F = [(1:6)',2*ones(6,1),-F_sam];
                Ele_pro_tmp(1:11,3) = A_sam(1);
                Ele_pro_tmp(12:23,3) = A_sam(2);
                Ele_pro_tmp(1:11,4) = E_sam(1);
                Ele_pro_tmp(12:23,4) = E_sam(2);

                fem = FEM(Node,Ele,Ele_pro_tmp,BC,F);
                u = fem.L_Slv_Sta;
                L(isam) = 0.14 - abs(u(20));
            end
        end

    end
end