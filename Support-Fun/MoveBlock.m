function [var_low,var_up] = MoveBlock(log_w)
% Moving-block variance estimator (Appendix B, Eqns. B.5-B.10):
%   estimates Var of a chain-wise mean estimator from Nc Markov chains.
%   var_est = R_eff * sigma_C^2 / Nc   (Eqn. B.8)
if iscell(log_w)
    Nite = length(log_w);
else
    Nite = 1;
end
log_mean = zeros(Nite,1);
var_est = zeros(Nite,1);
for ite = 1:Nite
    if Nite == 1
        log_mean(ite) = logmean(log_w,"all");
        w = exp(log_w);
    else
        log_mean(ite) = logmean(log_w{ite},"all");
        w = exp(log_w{ite});
    end
    C_hat = mean(w,2);                  % chain-wise sub-estimators C_hat_j (Eqn. B.5)
    Nc = size(C_hat,1);
    var_C = var(C_hat);                 % sigma_C^2 (Eqn. B.6)
    Reff = zeros(round(Nc/4)-1,1);
    for b = 2:round(Nc/4)
        N_block = Nc-b+1;
        S_block = zeros(N_block,1);     % moving block averages S_hat_j (Eqn. B.7)
        for i = 1:N_block
            S_block(i) = mean(C_hat(i:i+b-1));
        end
        var_S = var(S_block)*(N_block-1)/N_block;   % sigma_S^2 (Eqn. B.10)
        Reff(b-1) = b*var_S/var_C;                  % b*sigma_S^2/sigma_C^2 (Eqn. B.9)
    end
    Reff_max = max([Reff;1]);
    if log_mean(ite) == -inf
        var_est(ite) = 0;
    else
        var_est(ite) = Reff_max*var_C/Nc;           % Var(H_hat_m) = R_eff*sigma_C^2/Nc (Eqn. B.8)
    end
end
var_low = sum(var_est);
var_up = sum(sqrt(var_est)).^2;
end

