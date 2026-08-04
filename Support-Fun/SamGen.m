function [u,y,Ncal,baystc] = SamGen(baystc,N)
% generate i.i.d samples from q0
if isfield(baystc.G,'L')
    L = baystc.G.L;
else
    L = eye(baystc.Ndim);
end
u = L*randn(baystc.Ndim,N);
baystc = baystc.EvlY(u);
y = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
Ncal = baystc.Ncal;
end