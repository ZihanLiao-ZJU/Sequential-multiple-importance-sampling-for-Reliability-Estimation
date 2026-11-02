function b=configure(b)
defaults=struct('NumSam',1000,'MaxIte',10000,'Verbose',true, ...
 'stop_before_seeds',false,'select_at_cap',false,'level_offset',0,'schedule',{{}});
if isprop(b,'Nsam'),defaults.NumSam=b.Nsam;end
if ~isfield(b.G,'p'),b.G.p=0.1;end
for key=fieldnames(defaults)'
 name=key{1};if ~isfield(b.G,'run')||~isfield(b.G.run,name),b.G=semcs.setfield(b.G,['run.' name],defaults.(name));end
end
opts=semcs.sampler_defaults();
for key=fieldnames(opts)'
 name=key{1};if ~isfield(b.G,'sampler')||~isfield(b.G.sampler,name),b.G=semcs.setfield(b.G,['sampler.' name],opts.(name));end
end
end
