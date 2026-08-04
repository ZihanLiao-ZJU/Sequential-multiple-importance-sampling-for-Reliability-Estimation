function [var_low,var_up] = MoveBlock(estimate)
if iscell(estimate)
    Nite = length(estimate);
else
    Nite = 1;
end
z = zeros(Nite,1);
varpi = zeros(Nite,1);
for ite = 1:Nite
    if Nite == 1
        z(ite) = logmean(estimate,"all");
        pi = exp(estimate);
    else
        z(ite) = logmean(estimate{ite},"all");
        pi = exp(estimate{ite});
    end
    C = mean(pi,2);
    Nc = size(C,1);
    sigmaC = var(C);
    Reff = zeros(round(Nc/4)-1,1);
    for b = 2:round(Nc/4)
        N_M = Nc-b+1;
        B = zeros(N_M,1);
        for i = 1:N_M
            B(i) = mean(C(i:i+b-1));
        end
        sigmaB = var(B)*(N_M-1)/N_M;
        Reff(b-1) = b*sigmaB/sigmaC;
    end
    Reffmax = max([Reff;1]);
    if z(ite) == -inf
        varpi(ite) = 0;
    else
        varpi(ite) = Reffmax*sigmaC/Nc;
    end
end
var_low = sum(varpi);
var_up = sum(sqrt(varpi)).^2;
end

