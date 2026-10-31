classdef SuS
    % Scalar subset-simulation intermediate distribution (no iteration loop).
    % The two reliability toolboxes use this identical class.
    % CoordinateSystem: 'standard_normal' (PF.Dist) or 'physical' (PF.Rdis).
    % ThresholdRule: 'hazen' quantile or 'midpoint' at round(N*p).
    % ClipThreshold selects whether the final threshold is clamped at zero.
    properties
        TF
        Nfun
        Ndim
        Ncal = 0;
        X
        Y
        G = inf;
        H = nan;
        p = 0.1;
        Nsam = 1000;
        CoordinateSystem = 'standard_normal';
        ThresholdRule = 'hazen';
        ClipThreshold = true;
    end
    properties(Dependent)
        FlgCvg
    end
    methods
        function obj = SuS(PF,coordinate_system,threshold_rule,clip_threshold)
            obj.TF=PF;obj.Nfun=PF.Nfun+3;obj.Ndim=PF.Ndim;
            if nargin>=2,obj.CoordinateSystem=coordinate_system;end
            if nargin>=3,obj.ThresholdRule=threshold_rule;end
            if nargin>=4,obj.ClipThreshold=clip_threshold;end
            assert(any(strcmp(obj.CoordinateSystem,{'standard_normal','physical'})));
            assert(any(strcmp(obj.ThresholdRule,{'hazen','midpoint'})));
        end
        function obj = EvlY(obj,u)
            obj.X=u;
            if strcmp(obj.CoordinateSystem,'physical')
                f=obj.TF.EvlLSF(u);
                log_phi=obj.TF.Rdis.EvlPDF(u);
            else
                f=obj.TF.EvlLSF(obj.TF.Dist.U2X(u));
                log_phi=logGauss(u);
            end
            obj.Y=[log_phi;f];obj.Ncal=size(u,2);
        end
        function value = EvlLKF(obj)
            f=obj.Y(2,:);
            if strcmp(obj.CoordinateSystem,'physical') && isinf(obj.G)
                value=zeros(size(f));
            else
                value=log(double(f<obj.G));
            end
        end
        function value = EvlPDF(obj,varargin)
            if ~isempty(varargin)
                if strcmp(obj.CoordinateSystem,'physical')
                    value=obj.TF.Rdis.EvlPDF(obj.X,varargin{1});
                else
                    value=-0.5*(log(2*pi)+obj.X.^2);
                end
            elseif ~isempty(obj.Y)
                value=obj.Y(1,:);
            elseif strcmp(obj.CoordinateSystem,'physical')
                value=obj.TF.Rdis.EvlPDF(obj.X);
            else
                value=logGauss(obj.X);
            end
        end
        function [y,g,obj,N_cal] = UpdObj(obj,~,y,p)
            if nargin<4,p=obj.p;end
            f=y(4,:);
            if strcmp(obj.ThresholdRule,'midpoint')
                l=sort(f,'ascend');n=round(numel(f)*p);g=(l(n)+l(n+1))/2;
            else
                g=quantile(f,p);
            end
            if obj.ClipThreshold,g=max(g,0);end
            obj.G=g;N_cal=0;
        end
        function [x_seed,y_seed] = SltSed(obj,x,y,N)
            idx=y(4,:)<obj.G;
            n=sum(idx,'all');
            if n<=0,x_seed=[];y_seed=[];return;end
            ns=round(N./(1:N));nc=round(N./ns);
            nc=nc(find(n-nc>=0,1,'last'));
            if n>=nc,sel=randperm(n,nc);else,sel=randi(n,1,nc);end
            xx=x(:,idx);yy=y(:,idx);x_seed=xx(:,sel);y_seed=yy(:,sel);
            if ~isempty(y_seed)
                obj.X=x_seed;obj.Y=y_seed(3:end,:);
                y_seed=[obj.EvlLKF;obj.EvlPDF;obj.Y];
            end
        end
        function x = SamGen(obj,N)
            if strcmp(obj.CoordinateSystem,'physical')
                x=obj.TF.Rdis.SamGen(N);
            else
                x=randn(obj.Ndim,N);
            end
        end
        function value = get.FlgCvg(obj)
            value=obj.G<=0;
        end
    end
end
