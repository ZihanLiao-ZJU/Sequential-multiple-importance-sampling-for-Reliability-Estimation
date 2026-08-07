function log_Pi_hat = MisInt(log_alpha,log_W)
% Component-wise MIS estimator in log scale (Eqn. 10 / 33):
%   log_Pi_hat(ite) = log P_hat_ite = log[ (1/N) sum_k W(U_k^(ite)) ]
% where log_alpha = balance-heuristic weight (Eqn. 6) and
%       log_W     = log[ phi_n * I_F / q_hat ] part of the MIS weight W (Eqn. 8).
Nite = length(log_alpha);
log_Pi_hat = zeros(Nite,1);
for ite = 1:Nite
    log_Pi_hat(ite) = logmean(log_alpha{ite} + log_W{ite},"all");
end
end

