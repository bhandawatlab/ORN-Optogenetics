function dist = JSDiv(P,Q)
% JSDiv  Calculates the  jenson shannon divergence. 
%
% Inputs:
%    P: probability distribution
%    Q: probability distribution
%
% Outputs:
%    dist: JS divergence

M = 0.5*(P+Q);
dist = 0.5.*KLDiv(P,M)+0.5.*KLDiv(Q,M);

end