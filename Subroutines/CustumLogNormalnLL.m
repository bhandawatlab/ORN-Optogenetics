function [nLL] = CustumLogNormalnLL(y,pAll,pCovAll,nMu,nSyn,p)
% CustumLogNormalnLL  logNormal negative likelihood function
%
%   Inputs: y = vector of data points to evaluate
%           pAll = 
%           pCovAll = 
%           nMu = 
%           nSyn = 
%           p = 
%
%   Output: nLL = negative loglikelihood
%   

% 
% p_mu = (p(1:17));%11+6
% p_sig = (p(18:34));%11+6
% p_synergy = (p(35:50));%11+5
% p_covar = (p(51:66));%11+5

p_mu = (p(1:nMu));%11+6
p_sig = (p(nMu+1:2*nMu));%11+6
p_synergy = (p(2*nMu+1:2*nMu+nSyn));%11+5
p_covar = (p(2*nMu+nSyn+1:2*nMu+2*nSyn));%11+5

muF = (p_mu*pAll'+p_synergy*pCovAll')';
varF = (p_sig*pAll'+2.*p_covar*pCovAll')';

pdf = (1./(sqrt(varF).*sqrt(2*pi)) .* exp(-((y-muF).^2 ./ (2*(varF)))))+eps;%+eps;
nLL = -sum(log(pdf+eps));


end