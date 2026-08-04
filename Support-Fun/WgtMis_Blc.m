function w_mis = WgtMis_Blc(intBay,x,y,g,c)
% weight of mutiple importance sampling (balance heuristic)

beta = 1;
% initialization
Nite = length(y);
w_mis = cell(Nite,1);
Nsam = zeros(Nite,1);
Nsze = zeros(Nite,2);
for ite = 1:Nite
    Nsze(ite,:) = size(y{ite},2:3);
    Nsam(ite) = prod(Nsze(ite,:));
    w_mis{ite} = zeros(Nsze(ite,:));
end

for ite = 1:Nite
    pdf_crs = zeros(Nite,Nsam(ite));
    for ite1 = 1:Nite
        intBay.G = g{ite1};
        intBay.X = x{ite}(1:end,:);
        intBay.Y = y{ite}(3:end,:);
        pdf_crs(ite1,:) = intBay.EvlLKF+intBay.EvlPDF-c(ite1);
    end
    % balance heuristic
    w_mis_tmp = beta*log(Nsam(ite)) + beta*pdf_crs(ite,:) - logsum(beta*pdf_crs+beta*log(Nsam),1);
    % reshape
    w_mis{ite}(:) = w_mis_tmp(:);
end
end