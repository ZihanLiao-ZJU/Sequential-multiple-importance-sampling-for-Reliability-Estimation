classdef SuS
    % Scalar subset-simulation intermediate distribution (no iteration loop).
    % The two reliability toolboxes use this identical class.
    % Adapter: 'standard_normal' uses PF.Dist; 'physical' uses PF.Rdis.
    % Both sample/store Gaussian u; the adapter maps u to physical LSF inputs.
    % ThresholdRule: 'hazen' quantile or 'midpoint' at round(N*p).
    % ClipThreshold selects whether the final threshold is clamped at zero.
    properties
        TF
        Nfun
        Ndim
        Ncal = 0;
        X
        Y
        G = struct("parameters",inf,"p",0.1);
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
            obj.TF=PF;obj.Nfun=PF.Nfun+3;obj.Ndim=PF.Ndim;obj.G.L=eye(obj.Ndim);
            if nargin>=2,obj.CoordinateSystem=coordinate_system;end
            if nargin>=3,obj.ThresholdRule=threshold_rule;end
            if nargin>=4,obj.ClipThreshold=clip_threshold;end
            assert(any(strcmp(obj.CoordinateSystem,{'standard_normal','physical'})));
            assert(any(strcmp(obj.ThresholdRule,{'hazen','midpoint'})));
        end
        function obj = EvlY(obj,u)
            obj.X=u;
            if strcmp(obj.CoordinateSystem,'physical')
                f=obj.TF.EvlLSF(obj.TF.Rdis.U2X(u));
                log_phi=-.5*(obj.Ndim*log(2*pi)+sum(u.^2,1));
            else
                f=obj.TF.EvlLSF(obj.TF.Dist.U2X(u));
                log_phi=logGauss(u);
            end
            obj.Y=[log_phi;f];obj.Ncal=size(u,2);
        end
        function value = EvlLKF(obj)
            f=obj.Y(2,:);
            if strcmp(obj.CoordinateSystem,'physical') && isinf(obj.G.parameters)
                value=zeros(size(f));
            else
                value=log(double(f<obj.G.parameters));
            end
        end
        function value = EvlPDF(obj,varargin)
            if ~isempty(varargin)
                value=-.5*(log(2*pi)+obj.X(varargin{1},:).^2);
            elseif ~isempty(obj.Y)
                value=obj.Y(1,:);
            else
                value=-.5*(obj.Ndim*log(2*pi)+sum(obj.X.^2,1));
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
            obj.G.parameters=g;obj.G.p=p;g=obj.G;N_cal=0;
        end
        function [x_seed,y_seed] = SltSed(obj,x,y,N)
            idx=y(4,:)<obj.G.parameters;
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
            x=randn(obj.Ndim,N);
        end
        function value = get.FlgCvg(obj)
            value=obj.G.parameters<=0;
        end
    end
    methods

        function defaults = run_defaults(obj)
            defaults = struct('NumSam', 1000, 'ConPrb', 0.1, 'MaxIte', 10000, 'Verbose', true);
            if isprop(obj, 'Nsam')
                defaults.NumSam = obj.Nsam;
            end

            if isfield(obj.G, 'p')
                defaults.ConPrb = obj.G.p;
            end
        end

        function context = start_run(~, owner, sampler)
            context = SeMCS.start_context(owner, sampler, owner.MaxIte - 1, 'replace');
        end

        function context = advance_level(~, context, record, level)
            if level + 1 <= numel(context.owner.G)
                context.model.G = SeMCS.prescribed_state(context.model.G, context.owner.G{level+1}, false, 'replace');
                extra = 0;
            else
                [context.y{level}, ~, context.model, extra] = context.model.UpdObj(record.x, record.y, context.owner.ConPrb);
            end

            context.Ncal(level) = context.Ncal(level) + extra;
            context.g{level+1, 1} = context.model.G;
            context.sampler.BayStc = context.model;
        end

        function context = select_run_seeds(~, context, record, level)
            if context.targets(level + 1) <= 0
                error('SeMCS:ScheduleExhausted', 'Sample budget exhausted before seed selection.');
            end

            target = context.targets(level + 1);
            [context.xs, context.ys] = context.model.SltSed(record.x, context.y{level}, target);
            context.nc = size(context.xs, 2);
            context.steps = round(target / context.nc);
            if context.owner.Verbose
                fprintf('Level %d: N=%d, chains=%d, steps=%d, calls=%d\n', level, record.count, record.shape(1), record.shape(2), context.Ncal(level));
            end
        end

        function [yes, context] = run_complete(~, context, ~, phase)
            yes = false;
            if ~strcmp(phase, 'after_seeds')
                return;
            end

            if context.nc <= 0 || numel(context.x) == context.limit
                yes = true;
                return;
            end

            context.complete = context.model.FlgCvg;
            yes = context.complete;
        end

        function out = finish_run(~, context)
            out = SeMCS.pack_output(context);
        end

    end
end
