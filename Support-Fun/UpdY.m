function y = UpdY(IntBay,u,y)
% evaluate the likelihood function values
% ----------------------------------------------------------------------
% SYNTAX:
% y = EvlLKF(obj,theta)
% ----------------------------------------------------------------------
% INPUTS:
% obj   : class constructed
% u     : seed samples in standard normal space             [Ndim,Nsam]
% y     : function list of samples                         [Nfun+3,Nsam]
% ----------------------------------------------------------------------
% OUTPUTS:
% y     : output function values                           [Nfun+3,Nsam]
%         --intermediate likelihood function log L_i (ISF)      [1,Nsam]
%         --intermediate prior PDF log pi_i (rPDF)              [1,Nsam]
%         --standard normal log-PDF log phi_n                   [1,Nsam]
%         --PF value g                                      [Nfun,Nsam]
% ----------------------------------------------------------------------

% extraction the y except for log L_i and log pi_i
y = y(3:end,:);
IntBay.X = u;
IntBay.Y = y;
log_pi = IntBay.EvlPDF;
log_L = IntBay.EvlLKF;
y = [log_L;log_pi;y];
end
