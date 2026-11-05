function B = logstd(A,dim)
% Logarithmic version of sum
% Syntax:
% -----------------------------------------------------------------------------------
% B = logmean(A,dim)
% -----------------------------------------------------------------------------------
%
% Inputs:
% -----------------------------------------------------------------------------------
% A    : logarithmic matrix
% dim  : returns the std (in logarithm) of the elements along dimension dim
% -----------------------------------------------------------------------------------
%
% Outputs:
% -----------------------------------------------------------------------------------
% B    : std (in logarithm) array
% -----------------------------------------------------------------------------------
sz = size(A);
logmu = logmean(A,dim);
if strcmpi(dim,'all')
    B = (logsum(2*real(logminus(A,logmu)),dim) - log(prod(sz,dim)-1))/2;
else
    B = (logsum(2*real(logminus(A,logmu)),dim) - log(sz(dim)-1))/2;
end
return