function [r] = sampleKinematics_raw040521(self,state,t,f,df,ttBase,ttInh,stdRatio)
m = self.model.params{state};

% mapping f/df to kinematic distribution
mu = zeros(size(t));
[~,tSlice]=max(t<(m.KNN.tt+m.KNN.tt(2)./2),[],2);
uSlice =unique(tSlice);
for s = 1:numel(uSlice)
    slice = uSlice(s);
    mu(tSlice==slice) = m.g2{slice}(df(tSlice==slice),f(tSlice==slice));
end
mu(t<0) = m.baseline(1);% when t<0, it means before first entry

if sum(ttBase>=0 & ttBase<2)>0 && state ==7
    a = 1;
end
mu(ttBase>=0 & ttBase<2) = m.baselineEarly(1);
mu(ttBase>=2) = m.baselineLate(1);

sigma = zeros(size(t));
% sampling kinematics from that distribution
if strcmpi(m.fitFun, 'logN')
    for i = 1:unique(tSlice)
        sigma(tSlice==i) = m.g2_std{i}(df(tSlice==i),f(tSlice==i))./stdRatio;
    end
    sigma(ttBase>=0 & ttBase<2) = m.baselineEarly(2);
    sigma(ttBase>=2) = m.baselineLate(2);
    r = lognrnd(mu,sigma);
elseif strcmpi(m.fitFun,'norm')
    for i = 1:unique(tSlice)
        sigma(tSlice==i) = m.g2_std{i}(df(tSlice==i),f(tSlice==i))./stdRatio;
    end
    sigma(ttBase>=0 & ttBase<2) = m.baselineEarly(2);
    sigma(ttBase>=2) = m.baselineLate(2);
    r = normrnd(mu,sigma);
elseif strcmpi(m.fitFun, 'exp')
    r = exprnd(mu);
end

r(t<0) = m.baselineDat(ceil(numel(m.baselineDat).*rand(sum(t<0),1)));

if any(~isnan(ttInh)) && ~isempty(m.inhibitionKin.a{1})
    for s = 1:numel(uSlice)
        slice = uSlice(s);
        a(tSlice==slice) = m.inhibitionKin.a{slice}(ttInh(tSlice==slice));
        b(tSlice==slice) = m.inhibitionKin.b{slice}(ttInh(tSlice==slice));
    end
    
%     a = m.inhibitionKin.a(ttInh);
%     b = m.inhibitionKin.b(ttInh);
    if strcmpi(m.inhibitionKin.fitFun,'Lognormal') || strcmpi(m.inhibitionKin.fitFun,'logN')
        %--------------------log normal
        sampInh = lognrnd(a,b);% in mm/s
    elseif strcmpi(m.inhibitionKin.fitFun,'beta')
        %--------------------beta
        sampInh = betarnd(a,b).*m.inhibitionKin.maxVal;
    end
    r(~isnan(sampInh)) = sampInh(~isnan(sampInh));
end

qBound = quantile([m.baselineDat;m.allDat],0.99);
outOfBounds = r>qBound;%max([m.baselineDat;m.allDat]
r(outOfBounds) = max(qBound);
outOfBounds = r<min([m.baselineDat;m.allDat]);
r(outOfBounds) = min([m.baselineDat;m.allDat]);


end