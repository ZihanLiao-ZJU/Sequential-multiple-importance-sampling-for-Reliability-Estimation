function [u,y,N_cal,baystc] = SamGen(baystc,N)
% Generate i.i.d. samples from the initial ISD q_0 = phi_n(.; sigma^2 I)
%   (Eqn. 21 with i = 0: no truncation), via u = L_pi * randn.
if isfield(baystc.G,'L_pi')
    L_pi = baystc.G.L_pi;
else
    L_pi = eye(baystc.Ndim);
end
u = L_pi*randn(baystc.Ndim,N);
baystc = baystc.EvlY(u);
y = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];   % [log L; log pi; log phi_n; g]
N_cal = baystc.Ncal;
end
