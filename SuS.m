classdef SuS
    % Subset simulation with Nataf transfomration for reliability problems
    properties
        % Name       Description                                       Type              Size
        % -----------------------------------------------------------------------------------
        % LSF        Limit state function                              /                [1,1]
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
        G = inf;
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
        function obj = SuS(LSF)
            % Initializatipon
            % ----------------------------------------------------------------------
            % SYNTAX:
            % obj = BayLSF(Tfun)
            % ----------------------------------------------------------------------
            % INPUTS:
            % Tfun : Target function
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % obj  : constructed class
            % ----------------------------------------------------------------------
            obj.TF = LSF;
            obj.Nfun = LSF.Nfun+3;
            obj.Ndim = LSF.Ndim;
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
            %         --LSF value f                                         [1,Nsam]
            % ----------------------------------------------------------------------
            x = U2X(obj,u);
            f = obj.TF.EvlLSF(x);
            logp = logGauss(u);
            obj.X = u;
            obj.Y = [logp;f];
            obj.Ncal = size(u,2);
        end

        function Li = EvlLKF(obj)
            % Evaluate the Intermediate likelihood Li
            % ----------------------------------------------------------------------
            % SYNTAX:
            % Li = EvlLi(obj,f)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : class constructed                                        [1,1]
            % y     : function list of samples                              [2,Nsam]
            %         --parameter distribution PDF                          [1,Nsam]
            %         --LSF value                                           [1,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % Li    : intermediate likelihood function Li                   [1,Nsam]
            % ----------------------------------------------------------------------
            y = obj.Y;
            f = y(2,:); % LSF value
            g = obj.G;
            Li = log(double(f<g));
        end

        function Pi = EvlPDF(obj)
            % Evaluate the Intermediate prior PDF pi
            % ----------------------------------------------------------------------
            % SYNTAX:
            % Li = EvlLi(obj,f)
            % ----------------------------------------------------------------------
            % INPUTS:
            % obj   : class constructed                                        [1,1]
            % y     : function list of samples                              [2,Nsam]
            %         --parameter distribution PDF                          [1,Nsam]
            %         --LSF value                                           [1,Nsam]
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % Li    : intermediate likelihood function Li                   [1,Nsam]
            % ----------------------------------------------------------------------
            y = obj.Y;
            logp = y(1,:);
            Pi = logp;
        end

        function [y,g,obj,Ncal] = UpdObj(obj,~,y)
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
            %         --LSF value f                                         [1,Nsam]
            % p     : conditional probability in SuS
            % ----------------------------------------------------------------------
            % OUTPUTS:
            % g     : updated parameters for intermediate likelihood function
            % h     : updated parameters for intermediate prior distribution
            % obj   : class updated
            % ----------------------------------------------------------------------

            % update g
            f = y(4,:);
            g = max(quantile(f,obj.p),0);
            % update obj
            obj.G = g;
            Ncal = 0;
        end

        function FlgCvg = get.FlgCvg(obj)
            FlgCvg = obj.G<=0;
        end

        function Pf = EvlPf(~,in)
            Nite = in.Nite;
            y = in.y;
            g = in.g;

            % Normalize y
            % --------------------------------------------------------------------------
            c = ones(Nite+1,1);
            for ite = 1:Nite
                fth = g{ite+1};
                yi = y{ite}(4,:);
                w = yi<fth;
                c(ite+1) = c(ite)*mean(w,"all");
            end
            % --------------------------------------------------------------------------
            Pf = c(Nite+1);
        end
    end
end