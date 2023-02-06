function [f_orco] = calcKNN_Habituation_PermutationTest(f_orco,states)
% calcKNN_Habituation_PermutationTest  performs a permutation test to
%   compare KNN distributions across time. Each state instance retains it's
%   f/df location, but the time since first entry is shuffled. The KNN
%   space is then calculated based on this shuffled dataset. this is
%   repeated 100 times to generate a confidence bound and the empirical KNN
%   space is compared to this confidence bound.
%
%   Inputs: f_orco = fly object
%           states = index of states to look at (e.g. 1=sharp turn speed)
%
%   Output: f_orco = fly object with the permutation test results saved as
%               a field called HabituationBootstrap in the respective
%               states
%   

nIt = 100;

for s = 1:numel(states)
    m = f_orco.model.params{states(s)};
    alldSpk = m.alldSpk;
    allSpk = m.allSpk;
    allTime = m.allTimeFE;
    fitFun = m.fitFun;
    allDat = m.allDat;
    mu_empirical = m.val(:);
    nData = numel(allDat);
    
    modelParams = m.KNN;
    
    xGrid = f_orco.model.TurnBias.XX(:,1)';
    yGrid = f_orco.model.TurnBias.YY(1,:);
    zGrid = modelParams.tt;
    [XX,YY,ZZ] = ndgrid(xGrid,yGrid,zGrid);
    YY = YY(:);
    ZZ = ZZ(:);
    XX_reshape = ZZ(:);
    
    pHat_mu = zeros(numel(XX),nIt);
    for it = 1:nIt
        allTime_shuffle=allTime(randperm(nData));
        pHat = calcKNN(allDat,alldSpk,allSpk,allTime_shuffle,XX_reshape,YY,ZZ,modelParams,fitFun);
        pHat_mu(:,it) = pHat(:,1);
    end
    
    N = sum(~isnan(pHat_mu),2);
    y_sem = nanstd(pHat_mu,[],2)./sqrt(N);
    yCI95 = nan(numel(XX),2);
    goodNdx = find(N>30);
    CI95 = zeros(numel(goodNdx),2);
    for i = 1:numel(goodNdx)
        CI95(i,:) = tinv([0.025 0.975], N(goodNdx(i))-1);
    end
    yCI95(goodNdx,:) = y_sem(goodNdx).*CI95;
    mu_CI95 = yCI95+nanmean(pHat_mu,2);
    
    h = -(mu_empirical<mu_CI95(:,1))+(mu_empirical>mu_CI95(:,2));
    h(isnan(mu_empirical)) = nan;
    f_orco.model.params{states(s)}.HabituationBootstrap.h = reshape(h,size(XX));
    f_orco.model.params{states(s)}.HabituationBootstrap.CI = mu_CI95;
    f_orco.model.params{states(s)}.HabituationBootstrap.nIt = nIt;
end

end

function [pHat] = calcKNN(allDat,alldSpk,allSpk,allTime,XX,YY,ZZ,modelParams,fitFun)

ratio = modelParams.ratio;
K = modelParams.K;
Thresh = modelParams.Thresh;

if isempty(alldSpk)
    pHat = nan(numel(XX),2);
    maskedDat = [];
else
    X = [alldSpk, allSpk, allTime]./ratio;
    
    Y = [XX,YY,ZZ]./ratio;
    [idx, D] = knnsearch(X,Y,'K',K);
    mask = D<Thresh;
    
    maskedDat = allDat(idx);
    maskedDat(~mask) = nan;
    
    try
        minCount = 15;
        
        if strcmpi(fitFun, 'logN')
            pHat = nan(size(maskedDat,1),2);
            p = nan(size(maskedDat,1),2);
            for i = 1:size(maskedDat,1)%6
                if sum(mask(i,:))>minCount
                    dat = maskedDat(i,:);
                    dat(isnan(dat))= [];
                    pHat(i,:) = lognfit(dat+eps);
                end
            end
            
        elseif strcmpi(fitFun, 'exp')
            pHat = nan(size(maskedDat,1),1);
            for i = 1:size(maskedDat,1)
                if sum(mask(i,:))>minCount
                    dat = maskedDat(i,:);
                    dat(isnan(dat))= [];
                    pHat(i,:) = expfit(dat);
                end
            end
        elseif strcmpi(fitFun, 'norm')
            pHat = nan(size(maskedDat,1),2);
            p = nan(size(maskedDat,1),2);
            for i = 1:size(maskedDat,1)%6
                if sum(mask(i,:))>minCount
                    dat = maskedDat(i,:);
                    dat(isnan(dat))= [];
                    [pHat(i,1),pHat(i,2)] = normfit(dat);
                end
            end
        end
    catch
        disp('error with fitFun name, use logN, exp, or norm');
    end
    
end
end