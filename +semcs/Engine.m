classdef Engine
    % The only sequential sampling loop. Algorithm rules belong to the model.
    properties
        Model
        Sampler
    end
    methods
        function obj=Engine(model,sampler)
            obj.Model=semcs.configure(model);obj.Sampler=sampler;
        end
        function out=run(obj)
            b=obj.Model;s=obj.Sampler;r=b.G.run;
            validateattributes(r.NumSam,{'numeric'},{'vector','integer','positive','finite'});
            validateattributes(r.MaxIte,{'numeric'},{'scalar','integer','positive','finite'});
            cap=max(1,r.MaxIte-r.level_offset);N=repmat(r.NumSam(1),cap+1,1);
            if numel(r.NumSam)>1,N(:)=0;N(1:numel(r.NumSam))=r.NumSam(:);end
            preset=r.schedule;
            if ~isempty(preset),b.G=merge(b.G,preset{1});end
            x=cell(cap,1);y=x;g=cell(cap+1,1);final=x;
            counts=zeros(cap,1);shape=zeros(cap,2);calls=zeros(cap,1);
            converged=false;g{1}=b.G;
            for k=1:cap
                s.BayStc=b;g{k}=b.G;
                if k==1
                    if ismethod(b,'initial_samples'),[x{k},y{k},calls(k),b]=b.initial_samples(N(k));
                    else,[x{k},y{k},calls(k),b]=s.SamGen(N(k),1);end
                else
                    [x{k},y{k},calls(k),b]=s.SamGen(xs,ys,steps);
                end
                shape(k,:)=[size(y{k},2),size(y{k},3)];counts(k)=prod(shape(k,:));
                final{k}=b.G;
                if k+1<=numel(preset)
                    b.G=merge(b.G,preset{k+1});extra=0;
                else
                    [y{k},~,b,extra]=b.advance(x{k},y{k},b.G.p);
                end
                calls(k)=calls(k)+extra;g{k+1}=b.G;
                if r.stop_before_seeds && b.stop(),converged=true;break;end
                if k==cap && ~r.select_at_cap,break;end
                if N(k+1)<=0,error('semcs:Schedule','Sample budget exhausted.');end
                [xs,ys,steps]=b.select_seeds(x{k},y{k},N(k+1));
                if r.Verbose,fprintf('Level %d: N=%d, chains=%d, steps=%d, calls=%d\n',k,counts(k),shape(k,1),shape(k,2),calls(k));end
                if isempty(xs)||k==cap,break;end
                if b.stop(),converged=true;break;end
            end
            out=struct('x',{x(1:k)},'y',{y(1:k)},'g',{g(1:k+1)}, ...
                'Nsam',counts(1:k),'Nsze',shape(1:k,:),'Ncal',calls(1:k), ...
                'Nite',k,'flg_cvg',converged,'intBay',b,'sampled_final_G',{final(1:k)}, ...
                'engine','semcs.Engine','schema_version',1);
        end
    end
end
function next=merge(old,next)
if ~isstruct(next),tmp=old;tmp.parameters=next;next=tmp;end
for key={'run','sampler','lambda','HMC_time'}
 name=key{1};if isfield(old,name),next.(name)=old.(name);end
end
end
