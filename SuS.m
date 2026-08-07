classdef SuS
    % Subset simulation with Nataf transformation for reliability problems.
    % SuS is a special, restricted case of the SeMIS framework (Section 3.2):
    % no variance inflation (sigma = 1) and no surrogate truncation; the ISF is
    % simply L_i(u) = I[g(u) < l_i].
    properties
        % Name       Description                                       Type              Size
        % -----------------------------------------------------------------------------------
        % PF        Performance function                              /                [1,1]
        % Nfun       Length of function list                           double           [1,1]
        % Ndim       Number of dimension of parameters                 double           [1,1]
        % X          Samples in standard normal space                  double     [Ndim,Nsam]
        % Y          Y value of U                                      doubl         [4,Nsam]
        % -----------------------------------------------------------------------------------
        TF
        Nfun
        Ndim
        Ncal
        X
        Y
        G = inf;      % PF threshold l_i (interface field kept for the SeMIS framework)
        p = 0.1;
        Nsam = 1000;
    end

    properties
        % Properties to be updated:
        % Name       Description                                         Type            Size
        % -----------------------------------------------------------------------------------
        % G          Parameters determine intermediate Bay inf           double         [1,2]
        % -----------------------------------------------------------------------------------

    end

    properties(Dependent)
        % Dependent properties:
        % Name        Description                                        Type            Size
        % -----------------------------------------------------------------------------------
        % FlgCvg      flag of convergence                                logical        [1,1]
        % -----------------------------------------------------------------------------------
        FlgCvg
    end

    methods
        function obj = SuS(PF)
            % Initialization
            % ----------------------------------------------------------------------
            % SYNTAX:
            % obj = SuS(PF)
            % ----------------------------------------------------------------------
            % INPUTS:
            % Tfun : Target function
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % obj  : constructed class
            % ----------------------------------------------------------------------
            obj.TF = PF;
            obj.Nfun = PF.Nfun+3;
            obj.Ndim = PF.Ndim;
        end

        function obj = EvlY(obj,u)
            % evaluate the Y values
            % ----------------------------------------------------------------------
            % SYNTAX:
            % y = EvlLKF(obj,theta)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : class constructed
            % u     : Samples in standard normal space                   [Ndim,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % obj   : class with updated Y
            % y     : function list of samples                              [2,Nsam]
            %         --parameter distribution PDF P                        [1,Nsam]
            %         --PF value g                                         [1,Nsam]
            % ----------------------------------------------------------------------
            x = U2X(obj,u);
            g_val = obj.TF.EvlLSF(x);
            log_phi = logGauss(u);
            obj.X = u;
            obj.Y = [log_phi;g_val];
            obj.Ncal = size(u,2);
        end

        function log_L = EvlLKF(obj)
            % Evaluate the Intermediate likelihood Li = I[g(u) < l_i]
            % ----------------------------------------------------------------------
            % SYNTAX:
            % Li = EvlLi(obj,f)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : class constructed                                        [1,1]
            % y     : function list of samples                              [2,Nsam]
            %         --parameter distribution PDF                          [1,Nsam]
            %         --PF value                                           [1,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % Li    : intermediate likelihood function Li                   [1,Nsam]
            % ----------------------------------------------------------------------
            y = obj.Y;
            g_val = y(2,:);       % PF value g(u)
            l_th = obj.G;         % threshold l_i
            log_L = log(double(g_val < l_th));
        end

        function log_pi = EvlPDF(obj)
            % Evaluate the Intermediate prior PDF pi = phi_n(u) (no inflation, sigma = 1)
            % ----------------------------------------------------------------------
            % SYNTAX:
            % Li = EvlLi(obj,f)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : class constructed                                        [1,1]
            % y     : function list of samples                              [2,Nsam]
            %         --parameter distribution PDF                          [1,Nsam]
            %         --PF value                                           [1,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % Li    : intermediate likelihood function Li                   [1,Nsam]
            % ----------------------------------------------------------------------
            y = obj.Y;
            log_phi = y(1,:);
            log_pi = log_phi;
        end

        function [y,l,obj,N_cal] = UpdObj(obj,~,y)
            % Update obj with Y
            % ----------------------------------------------------------------------
            % SYNTAX:
            % [g,h,obj] = UpdObj(obj,y)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : class constructed                                        [1,1]
            % y     : output function values                                [4,Nsam]
            %         --intermediate likelihood function Li                 [1,Nsam]
            %         --intermediate prior PDF pi                           [1,Nsam]
            %         --parameter distribution PDF P                        [1,Nsam]
            %         --PF value g                                         [1,Nsam]
            % p     : conditional probability in SuS
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % l     : updated threshold for intermediate likelihood function
            % obj   : class updated
            % ----------------------------------------------------------------------

            % update the threshold l_i (Eqn. 24, SuS case)
            g_val = y(4,:);
            l = max(quantile(g_val,obj.p),0);
            % update obj
            obj.G = l;
            N_cal = 0;
        end

        function FlgCvg = get.FlgCvg(obj)
            FlgCvg = obj.G<=0;
        end

        function Pf_hat = EvlPf(~,in)
            Nite = in.Nite;
            y = in.y;
            l = in.g;              % thresholds l_i for each level

            % Normalize y: z_hat = prod_m H_hat_m (Eqn. 28 with H_hat_m of Eqn. 31)
            % --------------------------------------------------------------------------
            z_hat = ones(Nite+1,1);
            for ite = 1:Nite
                l_th = l{ite+1};
                yi = y{ite}(4,:);
                ind_H = yi < l_th;                       % I[g < l_{ite+1}] (Eqn. 31)
                z_hat(ite+1) = z_hat(ite)*mean(ind_H,"all");
            end
            % --------------------------------------------------------------------------
            Pf_hat = z_hat(Nite+1);
        end
    end
end
