function [u_seed,y_seed] = SedSlt(IntBay,u,y,N)
% Select the seeds for the next iteration (Eqn. 24):
%   S_{i+1} = {k : g(U_k^(i)) < l_{i+1} }  (implemented by weighted
%   resampling with weight proportional to the new ISF L_{i+1}).
% INPUTS:
%   IntBay : intermediate Bayesian model
%   u      : samples in standard normal space                     [Ndim,N]
%   y      : function list of u                                     [Nfun,N]
%   N      : sample size of the next iteration
% OUTPUTS:
%   u_seed : seed samples
%   y_seed : corresponding function list of u_seed
y = y(1:end,:);
u = u(1:end,:);
n_sam = size(y,2);
% weight for random resampling (proportional to the new ISF)
log_L_new = ResRat(IntBay,u,y,IntBay.G);
w_seed = exp(log_L_new);
Nc_eff = round(sum(w_seed,"all")^2 / sum(w_seed.^2,"all"));   % ESS
if Nc_eff<=0 || isnan(Nc_eff)
    u_seed = [];
    y_seed = [];
    return
else
    if all(w_seed(:)==0 | w_seed(:)==1)
        idx_can = logical(w_seed);
    else
        idx_can = randsample(n_sam,Nc_eff,true,w_seed);
    end
end
u_can = u(:,idx_can);
y_can = y(:,idx_can);
% numbers of required seeds
Ns = round(N./(1:N));
Nc_seed = round(N./Ns);
Nc_seed = Nc_seed(find(Nc_eff-Nc_seed>=0,1,"last"));
% generate seeds
if Nc_eff >= Nc_seed
    idx_seed = randperm(Nc_eff,Nc_seed);
else
    idx_seed = randi(Nc_eff,1,Nc_seed);
end
u_seed = u_can(:,idx_seed);
y_seed = y_can(:,idx_seed);
if ~isempty(y_seed)
    y_seed = UpdY(IntBay,u_seed,y_seed);
end
end
