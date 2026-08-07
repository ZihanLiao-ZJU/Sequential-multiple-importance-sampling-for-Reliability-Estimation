function [log_alpha, log_alpha_full] = WgtMis_Blc(intBay, u, y, g, log_z)
% Weight of multiple importance sampling (balance heuristic, Eqn. 6), in log scale.
%   log_alpha{ite}      : log alpha_ite(u) reshaped to [Nc,Ns]   (back-compat)
%   log_alpha_full{ite} : Nite x Nsam matrix; row j = log alpha_j(u), the
%                         balance-heuristic weight of the j-th ISD at every
%                         sample of level ite  (used by the (C.9) estimates)

beta = 1;
% initialization
Nite = length(y);
log_alpha = cell(Nite,1);
log_alpha_full = cell(Nite,1);
N = zeros(Nite,1);
Nc_Ns = zeros(Nite,2);
for ite = 1:Nite
    Nc_Ns(ite,:) = size(y{ite},2:3);
    N(ite) = prod(Nc_Ns(ite,:));
    log_alpha{ite} = zeros(Nc_Ns(ite,:));
end

for ite = 1:Nite
    log_q = zeros(Nite,N(ite));
    for ite1 = 1:Nite
        intBay.G = g{ite1};
        intBay.X = u{ite}(1:end,:);
        intBay.Y = y{ite}(3:end,:);
        log_q(ite1,:) = intBay.EvlLKF+intBay.EvlPDF-log_z(ite1);
    end
    % balance heuristic (Eqn. 6): full matrix, rows = ISDs
    log_alpha_mat = beta*log(N) + beta*log_q - logsum(beta*log_q+beta*log(N),1);
    log_alpha_full{ite} = log_alpha_mat;
    % reshape the ite-th row (back-compat with P_hat_i = logmean(log_alpha+log_W))
    log_alpha{ite}(:) = log_alpha_mat(ite,:);
end
end
