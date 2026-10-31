classdef SeMCS
    % Sequential Monte Carlo iteration, independent of the intermediate model.
    % Set sampler.BayStc to GSS_MIS, SuS or a SeMIS distribution object.
    properties
        Sampler
        IntBay
        NumSam = 1000;
        ConPrb = 0.1;
        MaxIte = 10000;
        G = {};
        H = {};
        Verbose = true;
    end
    methods
        function obj = SeMCS(sampler)
            obj.Sampler = sampler;
            obj.IntBay = sampler.BayStc;
            if isprop(obj.IntBay,'Nsam'), obj.NumSam = obj.IntBay.Nsam; end
            if isstruct(obj.IntBay.G) && isfield(obj.IntBay.G,'p')
                obj.ConPrb = obj.IntBay.G.p;
            elseif isprop(obj.IntBay,'p')
                obj.ConPrb = obj.IntBay.p;
            end
        end
        function out = RunIte(obj)
            intBay=obj.IntBay; sampler=obj.Sampler; maxi=obj.MaxIte;
            Nsam=zeros(maxi+1,1); Nsze=zeros(maxi+1,2); Ncal=zeros(maxi+1,1);
            x=cell(maxi+1,1);y=x;g=x;h=x;
            if numel(obj.NumSam)>1
                Nsam(1:numel(obj.NumSam))=obj.NumSam;
            else
                Nsam(:)=obj.NumSam;
            end
            Nsze(1,:)=[Nsam(1),1];
            g{1}=intBay.G;
            has_h=isprop(intBay,'H');
            if has_h, h{1}=intBay.H; end
            n_predef=numel(obj.G);
            if n_predef>0, g(1:n_predef)=obj.G; end
            if ~isempty(obj.H), h(1:numel(obj.H))=obj.H; end
            if isstruct(g{1}), intBay.G = g{1}; else, intBay.G.parameters=g{1};g{1}=intBay.G;end
            if has_h, intBay.H = h{1}; end
            sampler.BayStc = intBay;
            flg_cvg=false;
            for ite=1:maxi-1
                if ite==1
                    [x{ite},y{ite},Ncal(ite),intBay]=sampler.SamGen(Nsze(ite,1),Nsze(ite,2));
                else
                    [x{ite},y{ite},Ncal(ite),intBay]=sampler.SamGen(x_seed,y_seed,Nsze(ite,2));
                end
                Nsam(ite)=prod(Nsze(ite,:));
                dNcal=0;
                if ite+1<=n_predef
                    previous_state=intBay.G;
                    if isstruct(g{ite+1}),intBay.G=g{ite+1};else,intBay.G.parameters=g{ite+1};end
                    for key={'lambda','HMC_time'}
                        if isfield(previous_state,key{1}),intBay.G.(key{1})=previous_state.(key{1});end
                    end
                    g{ite+1}=intBay.G;
                    if has_h && ~isempty(h{ite+1}), intBay.H=h{ite+1}; end
                else
                    [y{ite},g{ite+1},intBay,dNcal]=intBay.UpdObj(x{ite},y{ite},obj.ConPrb);
                end
                if has_h, h{ite+1}=intBay.H; end
                Ncal(ite)=Ncal(ite)+dNcal;
                sampler.BayStc=intBay;
                if numel(obj.NumSam)>1 && ite+1>numel(obj.NumSam)
                    error('SeMCS:ScheduleExhausted', 'NumSam schedule exhausted before seed selection.');
                end
                if ismethod(intBay,'SltSed')
                    [x_seed,y_seed]=intBay.SltSed(x{ite},y{ite},Nsam(ite+1));
                else
                    [x_seed,y_seed]=SedSlt(intBay,x{ite},y{ite},Nsam(ite+1));
                end
                Nsze(ite+1,:)=[size(x_seed,2),round(Nsam(ite+1)/size(x_seed,2))];
                if obj.Verbose, DspRst(ite,Nsam(ite),Nsze(ite,:),Ncal(ite)); end
                if Nsze(ite+1,1)<=0 || ite+1==maxi
                    break
                elseif intBay.FlgCvg
                    flg_cvg=true;break
                end
            end
            out=struct('x',{x(1:ite)},'y',{y(1:ite)},'g',{g(1:ite+1)}, ...
                'Nsam',Nsam(1:ite),'Nsze',Nsze(1:ite,:),'Ncal',Ncal(1:ite), ...
                'Nite',ite,'flg_cvg',flg_cvg,'intBay',intBay);
            if has_h, out.h=h(1:ite+1); end
        end
    end
end
