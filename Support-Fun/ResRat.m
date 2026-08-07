function [log_L_new,log_max] = ResRat(intBay,u,y,g)
% Selection weight (log scale) of existing samples under the NEW ISF L_new:
%   log_L_new = log L_new(u) + log pi_new(u) - log L_old(u) - log pi_old(u)
% Since the rPDF pi is unchanged across iterations, this reduces to
%   log L_new(u) (samples lie in the support of L_old).
% reshape y
u = u(1:end,:);
y = y(1:end,:);
intBay.G = g;
intBay.X = u;
intBay.Y = y(3:end,:);
log_L = intBay.EvlLKF;
log_pi = intBay.EvlPDF;
log_L_new_tmp = log_L + log_pi - y(1,:) - y(2,:);
log_max = max(log_L_new_tmp,[],"all");
log_L_new_tmp = log_L_new_tmp - log_max;
log_L_new_tmp(isnan(log_L_new_tmp)) = -inf;
% reshape log_L_new
log_L_new = zeros(size(y,2:3));
log_L_new(:) = log_L_new_tmp;
end
