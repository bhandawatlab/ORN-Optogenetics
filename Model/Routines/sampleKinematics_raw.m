function [r] = sampleKinematics_raw(self,state,t,f,df,tt,tSinceInhibition,stdRatio)
m = self.model.params{state};

y = m.h(m.tt,m.params.y,t);
phi = m.h(m.tt,m.params.phi,t);
a1 = m.h(m.tt,m.params.a1,t);
a2 = m.h(m.tt,m.params.a2,t);
b1 = m.h(m.tt,m.params.b1,t);
b2 = m.h(m.tt,m.params.b2,t);

% mapping f/df to kinematic distribution
%mu = m.g(m.params.c1,m.params.c2,y,phi,a1,a2,b1,b2,f,df);
mu = zeros(size(t));
[~,tSlice]=max(t<m.tt+m.tt(2)./2,[],2);
for i = 1:unique(tSlice)
    mu(tSlice==i) = m.g2{i}(df(tSlice==i),f(tSlice==i));
end
mu(t<0) = m.baseline(1);% when t<0, it means before first entry

mu(tt>0 & tt<=2) = m.baselineEarly(1);
mu(tt>2) = m.baselineLate(1);

sigma = zeros(size(t));
% sampling kinematics from that distribution
if strcmpi(m.fitFun, 'logN')
%     sigma = m.sigma(mu,m.b);
%     sigma(t<0) = m.baseline(2);
%     sigma(tt>0 & tt<=2) = m.baselineEarly(2);
%     sigma(tt>2) = m.baselineLate(2);
    for i = 1:unique(tSlice)
        sigma(tSlice==i) = m.g2_std{i}(df(tSlice==i),f(tSlice==i))./stdRatio;
    end

    %sigma = m.g2_std{tSlice}(df,f)./stdRatio;
    %sigma = zeros(size(mu));
    r = lognrnd(mu,sigma);
elseif strcmpi(m.fitFun,'norm')
    for i = 1:unique(tSlice)
        sigma(tSlice==i) = m.g2_std{i}(df(tSlice==i),f(tSlice==i))./stdRatio;
    end
    %sigma = m.g2_std{tSlice}(df,f)./stdRatio;
    %sigma = zeros(size(mu));
    r = normrnd(mu,sigma);
elseif strcmpi(m.fitFun, 'exp')
    r = exprnd(mu);
end

% while any(outOfBounds)
%     r(outOfBounds) = lognrnd(mu(outOfBounds),sigma(outOfBounds));
%     outOfBounds = r>max(m.baselineDat);
% end

r(t<0) = m.baselineDat(ceil(numel(m.baselineDat).*rand(sum(t<0),1)));
%figure;histogram(r(t<0),[0.5:0.5:25])

if any(~isnan(tSinceInhibition))
    a = m.inhibitionKin.a(tSinceInhibition);
    b = m.inhibitionKin.b(tSinceInhibition);
    if strcmpi(m.inhibitionKin.fitFun,'Lognormal') || strcmpi(m.inhibitionKin.fitFun,'logN')
        %--------------------log normal
        sampInh = lognrnd(a,b);% in mm/s
    elseif strcmpi(m.inhibitionKin.fitFun,'beta')
        %--------------------beta
        sampInh = betarnd(a,b).*m.inhibitionKin.maxVal;
    end
    r(~isnan(sampInh)) = sampInh(~isnan(sampInh));
end

outOfBounds = r>max([m.baselineDat;m.allDat]);
r(outOfBounds) = max([m.baselineDat;m.allDat]);
outOfBounds = r<min([m.baselineDat;m.allDat]);
r(outOfBounds) = min([m.baselineDat;m.allDat]);


if any(isnan(r))
    a = 1;
end

end