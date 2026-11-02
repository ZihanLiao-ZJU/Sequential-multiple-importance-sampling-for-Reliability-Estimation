classdef MCMCBase
    properties(Dependent)
        SedSav
        Lambda
        AcpOpt
        NumIte
        Adapt
        Ntry
        InitialTime
    end

    % Common Gaussian/Nataf target protocol: prior times likelihood factor.
    properties
        BayStc
        
    end
    properties (Dependent)
        Ndim
        Nfun
    end
    methods
        function d=get.Ndim(obj), d=obj.BayStc.Ndim; end
        function n=get.Nfun(obj), n=obj.BayStc.Nfun; end
        function obj=MCMCBase(b)
            if nargin>0, obj.BayStc=b; end
        end
        function [b,u,y,n,save,L,Nc,Ns]=prepare(obj,varargin)
            b=semcs.configure(obj.BayStc); d=b.Ndim; L=eye(d);
            if isfield(b.G,'L'), L=b.G.L;
            elseif isfield(b.G,'L_pi'), L=b.G.L_pi; end
            if isscalar(L), L=L*eye(d); end
            if ~isequal(size(L),[d,d]) || any(~isfinite(L(:))) || rank(L)~=d
                error('MCMC:Reference','Gaussian factor must be finite nonsingular Ndim by Ndim.');
            end
            if numel(varargin)==2
                Nc=varargin{1}; Ns=varargin{2};
                if ~isscalar(Nc) || ~isscalar(Ns) || ~isfinite(Nc) || ~isfinite(Ns) || Nc<1 || Ns<1 || Nc~=floor(Nc) || Ns~=floor(Ns), error('MCMC:Count','Counts must be positive integers.'); end
                if ismethod(b,'SamGen'), u=b.SamGen(Nc); else, u=L*randn(d,Nc); end
                [b,y,n]=obj.evaluate(b,u); save=true;
            elseif numel(varargin)==3
                u=varargin{1}; y=varargin{2}; Ns=varargin{3}; Nc=size(u,2); n=0; save=b.G.sampler.SedSav;
            else
                error('MCMC:Arguments','SamGen(Nchains,Nsteps) or SamGen(U_seed,Y_seed,Nsteps).');
            end
            if Nc<1 || Ns<1 || Nc~=floor(Nc) || Ns~=floor(Ns) || ~isequal(size(u),[d,Nc]) || ~isequal(size(y),[b.Nfun,Nc])
                error('MCMC:Shape','Invalid chain count, step count or packed seed outputs.');
            end
        end
        function [b,y,n]=evaluate(~,b,u)
            b=b.EvlY(u); y=[b.EvlLKF;b.EvlPDF;b.Y]; n=b.Ncal;
            if ~isequal(size(y),[b.Nfun,size(u,2)])
                error('MCMC:Output','Invalid packed target outputs.');
            end
        end
        function r=reject(~,new,old,p)
            delta=new-old; r=isnan(delta) | delta<log(max(p,realmin));
        end
        function value=state(~,b,key,default)
            value=default; if isfield(b.G,key) && ~isempty(b.G.(key)), value=b.G.(key); end
        end
    end
    methods
        function value=option(~,b,name),value=b.G.sampler.(name);end

    end

    methods
        function value=get.SedSav(obj)
            defaults=semcs.sampler_defaults();value=defaults.SedSav;
            if ~isempty(obj.BayStc),value=semcs.getfield(obj.BayStc.G,'sampler.SedSav',value);end
        end
        function obj=set.SedSav(obj,value)
            assert(~isempty(obj.BayStc),'semcs:BindFirst','Bind BayStc before setting sampler options.');
            obj.BayStc.G=semcs.setfield(obj.BayStc.G,'sampler.SedSav',value);
        end
        function value=get.Lambda(obj)
            defaults=semcs.sampler_defaults();value=defaults.Lambda;
            if ~isempty(obj.BayStc),value=semcs.getfield(obj.BayStc.G,'sampler.Lambda',value);end
        end
        function obj=set.Lambda(obj,value)
            assert(~isempty(obj.BayStc),'semcs:BindFirst','Bind BayStc before setting sampler options.');
            obj.BayStc.G=semcs.setfield(obj.BayStc.G,'sampler.Lambda',value);
        end
        function value=get.AcpOpt(obj)
            defaults=semcs.sampler_defaults();value=defaults.AcpOpt;
            if ~isempty(obj.BayStc),value=semcs.getfield(obj.BayStc.G,'sampler.AcpOpt',value);end
        end
        function obj=set.AcpOpt(obj,value)
            assert(~isempty(obj.BayStc),'semcs:BindFirst','Bind BayStc before setting sampler options.');
            obj.BayStc.G=semcs.setfield(obj.BayStc.G,'sampler.AcpOpt',value);
        end
        function value=get.NumIte(obj)
            defaults=semcs.sampler_defaults();value=defaults.NumIte;
            if ~isempty(obj.BayStc),value=semcs.getfield(obj.BayStc.G,'sampler.NumIte',value);end
        end
        function obj=set.NumIte(obj,value)
            assert(~isempty(obj.BayStc),'semcs:BindFirst','Bind BayStc before setting sampler options.');
            obj.BayStc.G=semcs.setfield(obj.BayStc.G,'sampler.NumIte',value);
        end
        function value=get.Adapt(obj)
            defaults=semcs.sampler_defaults();value=defaults.Adapt;
            if ~isempty(obj.BayStc),value=semcs.getfield(obj.BayStc.G,'sampler.Adapt',value);end
        end
        function obj=set.Adapt(obj,value)
            assert(~isempty(obj.BayStc),'semcs:BindFirst','Bind BayStc before setting sampler options.');
            obj.BayStc.G=semcs.setfield(obj.BayStc.G,'sampler.Adapt',value);
        end
        function value=get.Ntry(obj)
            defaults=semcs.sampler_defaults();value=defaults.Ntry;
            if ~isempty(obj.BayStc),value=semcs.getfield(obj.BayStc.G,'sampler.Ntry',value);end
        end
        function obj=set.Ntry(obj,value)
            assert(~isempty(obj.BayStc),'semcs:BindFirst','Bind BayStc before setting sampler options.');
            obj.BayStc.G=semcs.setfield(obj.BayStc.G,'sampler.Ntry',value);
        end
        function value=get.InitialTime(obj)
            defaults=semcs.sampler_defaults();value=defaults.InitialTime;
            if ~isempty(obj.BayStc),value=semcs.getfield(obj.BayStc.G,'sampler.InitialTime',value);end
        end
        function obj=set.InitialTime(obj,value)
            assert(~isempty(obj.BayStc),'semcs:BindFirst','Bind BayStc before setting sampler options.');
            obj.BayStc.G=semcs.setfield(obj.BayStc.G,'sampler.InitialTime',value);
        end

    end

end
