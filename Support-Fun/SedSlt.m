function [x_sed,y_sed] = SedSlt(IntBay,x,y,Nsam)
% select the seeds for next iteration
% Syntax:
% ----------------------------------------------------------------------
% [x_sed,y_sed] = SedSlt(IntL,IntP,x,y,p)
% ----------------------------------------------------------------------
% INPUTS:
% IntL  : Intermediate likelihood function
% IntP  : Intermediate prior distribution
% x     : generated samples from sampler                     [Ndim,Nsam]
% y     : target function value of u                            [6,Nsam]
%        --intermediate likelihood function Li
%        --likelihood function L
%        --target function f
%        --intermediate prior PDF pi
%        --prior PDF p
%        --reference distribution PDF q
% p     : filtering ratio
% ----------------------------------------------------------------------
% OUTPUTS:
% x_sed : seed samples
% y_sed : corresponding function list of x_sed
% ----------------------------------------------------------------------
y = y(1:end,:);
x = x(1:end,:);
N = size(y,2);
% weight for random resampling
lograt = ResRat(IntBay,x,y,IntBay.G);
w_rat = exp(lograt);
Nc_max = round(sum(w_rat,"all")^2 / sum(w_rat.^2,"all"));
if Nc_max<=0 || isnan(Nc_max)
    x_sed = [];
    y_sed = [];
    return
else
    if all(w_rat(:)==0 | w_rat(:)==1)
        ind_can = logical(w_rat);
    else
        ind_can = randsample(N,Nc_max,true,w_rat);
    end
end
x_can = x(:,ind_can);
y_can = y(:,ind_can);
% numbers of required seeds
Ns = round(Nsam./(1:Nsam));
Nc_can = round(Nsam./Ns);
Nc_can = Nc_can(find(Nc_max-Nc_can>=0,1,"last"));
% generate seeds
if Nc_max >= Nc_can
    idx = randperm(Nc_max,Nc_can);
else
    idx = randi(Nc_max,1,Nc_can);
end
x_sed = x_can(:,idx);
y_sed = y_can(:,idx);
if ~isempty(y_sed)
    y_sed = UpdY(IntBay,x_sed,y_sed);
end
end