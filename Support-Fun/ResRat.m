function [lograt,c] = ResRat(intBay,x,y,g)
% reshape y
y = y(1:end,:);
x = x(1:end,:);
intBay.G = g;
intBay.X = x;
intBay.Y = y(3:end,:);
Li = intBay.EvlLKF;
pi = intBay.EvlPDF;
lograt_tmp = Li+pi-y(1,:)-y(2,:);
c = max(lograt_tmp,[],"all");
lograt_tmp = lograt_tmp-c;
lograt_tmp(isnan(lograt_tmp)) = -inf;
% reshape lograt
lograt = zeros(size(y,2:3));
lograt(:) = lograt_tmp;
end