classdef SeMIS
    % =======================================================================================
    % Definition
    % ---------------------------------------------------------------------------------------
    % Sequential Multiple Importance Sampling (SeMIS) framework.
    %
    % The algorithm sequentially constructs a series of intermediate importance sampling
    % distributions (ISDs) q_i(u) = pi_i(u) * L_i(u) / z_i  (Eqn. 11/21) and generates
    % samples from each ISD to progressively approach the failure domain. At each iteration:
    %
    %   1) Samples are generated from the current ISD q_i.
    %   2) The intermediate distribution parameters (l_i, l_hat_i, ...) are updated.
    %   3) Seed samples S_{i+1} are selected for the next iteration (Eqn. 24).
    %   4) New samples are produced from seeds using pMCMC (Eqn. 22-23).
    %
    % The resulting sequence of ISDs and samples is used for multiple importance sampling
    % (MIS) estimation of the failure probability (Eqn. 9-10).
    %
    % ---------------------------------------------------------------------------------------
    % Configuration
    % ---------------------------------------------------------------------------------------
    % Developer:
    %   Zihan Liao, Zhejiang University
    %   Weili Xia,  Zhejiang University
    %   Xiao  He,   Zhejiang University
    %
    % Supervisor:
    %   Binbin Li,  Zhejiang University
    %
    % Initial Version : Jan 2023
    % Last Modified   : Jan 2026
    %
    % ---------------------------------------------------------------------------------------
    % Reference
    % ---------------------------------------------------------------------------------------
    % Sequential Multiple Importance Sampling framework for reliability estimation.
    % =======================================================================================

    properties (Constant)
        % properties:
        % Name      Description                         Type        Size
        % -----------------------------------------------------------------------
        % MaxIte    Maximum number of iterations        double      [1,1]
        MaxIte = 10000;
    end

    properties
        % properties:
        % Name        Description                                      Type
        % ---------------------------------------------------------------------------------
        % IntBay      Intermediate Bayesian model defining ISD         object
        % Sampler     Sampling object (must contain BayStc, SamGen)    object
        % NumSam      Number of samples per iteration (N in paper)     double / vector
        % G           Predefined parameters of intermediate ISDs       cell / vector
        % X           Storage for generated samples                    cell
        % Y           Storage for function evaluations                 cell

        % Intermediate Bayesian model
        IntBay
        
        % Sampler object (must contain .BayStc, .SamGen methods)
        Sampler

        % Number of samples per iteration N (scalar or vector)
        NumSam

        % Predefined intermediate distribution parameters
        G = [];

        % Storage for all samples
        X = [];
        Y = [];
    end

    methods
        function obj = SeMIS(IntBay,Sampler)
            % Constructor
            % ----------------------------------------------------------------------
            % SYNTAX:
            % obj = SeMIS(IntBay,Sampler)
            %
            % INPUTS:
            % IntBay   : intermediate Bayesian model defining ISDs
            % Sampler  : sampling object supporting BayStc and SamGen
            %
            % OUTPUTS:
            % obj      : constructed SeMIS object
            % ----------------------------------------------------------------------

            if nargin < 1
                error('Inputs must include the Sampler object with a valid BayStc field.');
            end
            obj.Sampler = Sampler;
            obj.IntBay = IntBay;
            obj.NumSam = IntBay.Nsam;
        end

        function out = RunIte(obj)
            % Perform sequential sampling iterations of SeMIS
            % ----------------------------------------------------------------------
            % SYNTAX:
            % out = RunIte(obj)
            %
            % DESCRIPTION:
            % Executes the iterative SeMIS procedure (Algorithm 1):
            %   1) Generate samples from the current ISD q_i.
            %   2) Update the intermediate distribution parameters.
            %   3) Select seed samples S_{i+1} for the next iteration.
            %   4) Generate samples from seeds using pMCMC.
            %
            % Iterations terminate when the convergence criterion defined in
            % the intermediate Bayesian model is satisfied (l_{i+1} <= 0).
            %
            % OUTPUTS:
            % out : structure containing
            %       x        - samples (standard-normal space U) from each iteration
            %       y        - function evaluations
            %       g        - intermediate distribution parameters
            %       Nsam     - number of samples per iteration (N)
            %       Nsze     - chain structure [#chains Nc, samples per chain Ns]
            %       Ncal     - number of function evaluations
            %       Nite     - number of iterations performed
            %       flg_cvg  - convergence flag
            %       intBay   - final intermediate Bayesian model
            % ----------------------------------------------------------------------

            % Extract settings
            intBay  = obj.IntBay;
            sampler = obj.Sampler;
            max_ite = obj.MaxIte;
            flg_N = numel(obj.NumSam) > 1;
            flg_predef = ~isempty(obj.G);

            % Initialize storage
            N = zeros(max_ite,1);        % sample size per iteration
            Nc_Ns = zeros(max_ite,2);    % chain structure [#chains, samples/chain]
            N_cal = zeros(max_ite,1);    % number of function evaluations
            u = cell(max_ite,1); y = cell(max_ite,1); g = cell(max_ite,1);

            % Sample size setup
            if flg_N
                n_predef = length(obj.NumSam);
                N(1:n_predef) = obj.NumSam;
            else
                N(:) = obj.NumSam;
            end
            Nc_Ns(1,:) = [N(1), 1];

            % Intermediate distribution setup
            if flg_predef
                n_predef = length(obj.G);
                g(1:n_predef) = obj.G;
            else
                n_predef = 0;
                g{1} = intBay.G;
            end

            % ===================== Sequential iterations =====================
            for ite = 1:max_ite
                % === Sample Generation ===
                if all(Nc_Ns(ite,:) > 0) && ~isnan(Nc_Ns(ite,2))
                    if ite == 1
                        [u{ite}, y{ite}, N_cal(ite), intBay] = SamGen(intBay, Nc_Ns(ite,1));
                    else
                        [u{ite}, y{ite}, N_cal(ite), intBay] = sampler.SamGen(u_seed, y_seed, Nc_Ns(ite,2));
                    end
                end
                N(ite) = prod(Nc_Ns(ite,:));

                % === Update Intermediate Distribution ===
                if ite + 1 <= n_predef
                    intBay.G = g{ite+1};
                else
                    [y{ite}, g{ite+1}, intBay, dN_cal] = intBay.UpdObj(u{ite}, y{ite});
                end
                sampler.BayStc = intBay;
                N_cal(ite) = N_cal(ite) + dN_cal;

                % === Seed Selection (Eqn. 24) ===
                % Select samples consistent with the next ISD
                [u_seed, y_seed] = SedSlt(intBay, u{ite}, y{ite}, N(ite+1));

                Nc_Ns(ite+1,1) = size(u_seed,2);                    % number of chains Nc
                Nc_Ns(ite+1,2) = round(N(ite+1) / Nc_Ns(ite+1,1));  % samples per chain Ns

                % === Display Iteration Summary ===
                DspRst(ite, N(ite), Nc_Ns(ite,:), N_cal(ite));

                % === Convergence Check ===
                if Nc_Ns(ite+1,1) <= 0 || ite + 1 == max_ite
                    flg_cvg = false;
                    break
                elseif intBay.FlgCvg
                    flg_cvg = true;
                    break
                end
            end

            % ===================== Post-processing =====================
            Nite = ite;
            [u, y, N, Nc_Ns, N_cal] = VabSrk(Nite, u, y, N, Nc_Ns, N_cal);
            g = VabSrk(Nite + 1, g);

            % ===================== Output =====================
            out.x = u;             % samples in the standard-normal space (paper: U)
            out.y = y;
            out.g = g;
            out.Nsam = N;
            out.Nsze = Nc_Ns;
            out.Ncal = N_cal;
            out.Nite = Nite;
            out.flg_cvg = flg_cvg;
            out.intBay = intBay;
        end
    end
end
