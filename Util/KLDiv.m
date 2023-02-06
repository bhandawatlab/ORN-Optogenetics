function dist = KLDiv(P,Q)
% KLDiv  Calculates the  Kullback-Leibler divergence
%
% Inputs:
%    P: probability distribution
%    Q: probability distribution
%
% Outputs:
%    dist: JS divergence

P = (P+eps)./sum(P+eps);
Q = (Q+eps)./sum(Q+eps);
dist = sum(P .* (log2(P) - log2(Q)), 1);

end