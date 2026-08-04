function y = UpdY(IntBay,x,y)
% evaluate the likelihood function values
% ----------------------------------------------------------------------
% SYNTAX:
% y = EvlLKF(obj,theta)
% ----------------------------------------------------------------------
% INPUTS:
% obj   : class constructed
% y     : function list of samples                         [Nfun+3,Nsam]
% ----------------------------------------------------------------------
% OUTPUTS:
% y     : output function values                           [Nfun+3,Nsam]
%         --intermediate likelihood function Li                 [1,Nsam]
%         --intermediate prior PDF pi                           [1,Nsam]
%         --parameter distribution PDF P                        [1,Nsam]
%         --LSF value f                                      [Nfun,Nsam]
% ----------------------------------------------------------------------

% extraction the y except for Li and pi
y = y(3:end,:);
IntBay.X = x;
IntBay.Y = y;
pi = IntBay.EvlPDF;
Li = IntBay.EvlLKF;
y = [Li;pi;y];
end