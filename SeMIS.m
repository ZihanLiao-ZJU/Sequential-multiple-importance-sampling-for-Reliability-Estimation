classdef SeMIS
    % =======================================================================================
    % Definition
    % ---------------------------------------------------------------------------------------
    % Sequential Multiple Importance Sampling (SeMIS) framework.
    %
    % The algorithm sequentially constructs a series of intermediate importance sampling
    % distributions (ISDs) and generates samples from each ISD to progressively approach
    % the failure domain. At each iteration:
    %
    %   1) Samples are generated from the current ISD.
    %   2) The intermediate distribution parameter is updated.
    %   3) Seed samples are selected for the next iteration.
    %   4) New samples are produced from seeds using MCMC.
    %
    % The resulting sequence of ISDs and samples is used for multiple importance sampling
    % (MIS) estimation of the failure probability.
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
        % NumSam      Number of samples per iteration                  double / vector
        % G           Predefined parameters of intermediate ISDs       cell / vector
        % X           Storage for generated samples                    cell
        % Y           Storage for function evaluations                 cell

        % Intermediate Bayesian model
        IntBay
        
        % Sampler object (must contain .BayStc, .SamGen methods)
        Sampler

        % Number of samples per iteration (scalar or vector)
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
            % Executes the iterative SeMIS procedure:
            %   1) Generate samples from the current ISD.
            %   2) Update the intermediate distribution parameter.
            %   3) Select seed samples for the next iteration.
            %   4) Generate samples from seeds using MCMC.
            %
            % Iterations terminate when the convergence criterion defined in
            % the intermediate Bayesian model is satisfied.
            %
            % OUTPUTS:
            % out : structure containing
            %       x        - samples from each iteration
            %       y        - function evaluations
            %       g        - intermediate distribution parameters
            %       Nsam     - number of samples per iteration
            %       Nsze     - chain structure [#chains, samples per chain]
            %       Ncal     - number of function evaluations
            %       Nite     - number of iterations performed
            %       flg_cvg  - convergence flag
            %       intBay   - final intermediate Bayesian model
            % ----------------------------------------------------------------------

            % Extract settings
            intBay  = obj.IntBay;
            sampler = obj.Sampler;
            maxi    = obj.MaxIte;
            flg_Nsam = numel(obj.NumSam) > 1;
            flg_IntD = ~isempty(obj.G);

            % Initialize storage
            Nsam = zeros(maxi,1);      % sample size per iteration
            Nsze = zeros(maxi,2);      % chain structure [#chains, samples/chain]
            Ncal = zeros(maxi,1);      % number of function evaluations
            x = cell(maxi,1); y = cell(maxi,1); g = cell(maxi,1);

            % Sample size setup
            if flg_Nsam
                Nite_sam = length(obj.NumSam);
                Nsam(1:Nite_sam) = obj.NumSam;
            else
                Nsam(:) = obj.NumSam;
            end
            Nsze(1,:) = [Nsam(1), 1];

            % Intermediate distribution setup
            if flg_IntD
                Nite_int = length(obj.G);
                g(1:Nite_int) = obj.G;
            else
                Nite_int = 0;
                g{1} = intBay.G;
            end

            % ===================== Sequential iterations =====================
            for ite = 1:maxi
                % === Sample Generation ===
                if all(Nsze(ite,:) > 0) && ~isnan(Nsze(ite,2))
                    if ite == 1
                        [x{ite}, y{ite}, Ncal(ite),intBay] = SamGen(intBay,Nsze(ite,1));
                    else
                        [x{ite}, y{ite}, Ncal(ite),intBay] = sampler.SamGen(x_sed, y_sed, Nsze(ite,2));
                    end
                end
                Nsam(ite) = prod(Nsze(ite,:));

                % === Update Intermediate Distribution ===
                if ite + 1 <= Nite_int
                    intBay.G = g{ite+1};
                else
                    [y{ite}, g{ite+1}, intBay, dNcal] = intBay.UpdObj(x{ite},y{ite});
                end
                sampler.BayStc = intBay;
                Ncal(ite) = Ncal(ite)+dNcal;

                % === Seed Selection ===
                % Select samples consistent with the next ISD
                [x_sed, y_sed] = SedSlt(intBay, x{ite}, y{ite}, Nsam(ite+1));

                Nsze(ite+1,1) = size(x_sed,2);                    % number of chains
                Nsze(ite+1,2) = round(Nsam(ite+1) / Nsze(ite+1,1)); % samples per chain

                % === Display Iteration Summary ===
                DspRst(ite, Nsam(ite), Nsze(ite,:), Ncal(ite));

                % === Convergence Check ===
                if Nsze(ite+1,1) <= 0 || ite + 1 == maxi
                    flg_cvg = false;
                    break
                elseif intBay.FlgCvg
                    flg_cvg = true;
                    break
                end
            end

            % ===================== Post-processing =====================
            Nite = ite;
            [x, y, Nsam, Nsze, Ncal] = VabSrk(Nite, x, y, Nsam, Nsze, Ncal);
            g = VabSrk(Nite + 1, g);

            % ===================== Output =====================
            out.x = x;
            out.y = y;
            out.g = g;
            out.Nsam = Nsam;
            out.Nsze = Nsze;
            out.Ncal = Ncal;
            out.Nite = Nite;
            out.flg_cvg = flg_cvg;
            out.intBay = intBay;
        end
    end
end