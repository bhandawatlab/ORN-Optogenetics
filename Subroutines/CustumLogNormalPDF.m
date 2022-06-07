function [pdf] = CustumLogNormalPDF(y,pAll,pCovAll,nMu,nSyn,p)

p_mu = cell2mat(p(1:nMu));%11+6
p_sig = cell2mat(p(nMu+1:2*nMu));%11+6
p_synergy = cell2mat(p(2*nMu+1:2*nMu+nSyn));%11+5
p_covar = cell2mat(p(2*nMu+nSyn+1:2*nMu+2*nSyn));%11+5


muF = (p_mu*pAll'+p_synergy*pCovAll')';
varF = (p_sig*pAll'+2.*p_covar*pCovAll')';

pdf = (1./(sqrt(varF).*sqrt(2*pi)) .* exp(-((y-muF).^2 ./ (2*(varF)))))+eps;%+eps;


% pdf = @(y,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22) ...
%     (1./(sqrt(varF([p7,p8,p9,p10,p11,p12],[p18,p19,p20,p21,p22])).*sqrt(2*pi)) .* ...
%     exp(-((y-muF([p1,p2,p3,p4,p5,p6],[p13,p14,p15,p16,p17])).^2 ./ (2*(varF([p7,p8,p9,p10,p11,p12],[p18,p19,p20,p21,p22]))))))+eps;%+eps;

%nLL = @(y,p) -sum(log(pdf(y,p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10),p(11),p(12),p(13),p(14),p(15),p(16),p(17),p(18),p(19),p(20),p(21),p(22))+eps));


end