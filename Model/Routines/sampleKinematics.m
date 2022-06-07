function [r] = sampleKinematics(self,state,t,f,df,tt)
m = self.model.params{state};

y = m.h(m.tt,m.params.y,t);
phi = m.h(m.tt,m.params.phi,t);
a1 = m.h(m.tt,m.params.a1,t);
a2 = m.h(m.tt,m.params.a2,t);
b1 = m.h(m.tt,m.params.b1,t);
b2 = m.h(m.tt,m.params.b2,t);

% mapping f/df to kinematic distribution
%mu = m.g(m.params.c1,m.params.c2,y,phi,a1,a2,b1,b2,f,df);
mu = m.g2(df,f);
mu(t<0) = m.baseline(1);% when t<0, it means before first entry

mu(tt>0 & tt<=2) = m.baselineEarly(1);
mu(tt>2) = m.baselineLate(1);

% sampling kinematics from that distribution
if strcmpi(m.fitFun, 'logN')
    %sigma = m.sigma(mu,m.b);
%     sigma = m.g2_std(df,f);
%     sigma(t<0) = m.baseline(2);
%     sigma(tt>0 & tt<=2) = m.baselineEarly(2);
%     sigma(tt>2) = m.baselineLate(2);
    
    sigma = zeros(size(mu));
    r = lognrnd(mu,sigma);
elseif strcmpi(m.fitFun, 'exp')
    r = exprnd(mu);
end


end