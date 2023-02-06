function [] = rulesOfSummationAnalysis(genAll,meta)
close all
tic;

t1 = -20;
t2 = 20;
t3 = 15;

%%
% % precompute bounds
% [stateKinAll,stateKinBaseline,stateKinBaselineMu,stateKinBaselineVar,...
%     StateLabel,f_orco] = getKinematicsByRegions(genAll,meta,t1,t2,t3,false);
% RegionLabel = {'NegDF','PosDF','HighF','Inh','LowF'};
% 
% %% Orco Ir8a
% lab = [{'Orco','Ir8a','OrcoIr8a','Orco','Ir8a'};
%  {'Before','Before','Before','',''}];
% 
% labCov = [{'Orco'};
%     {'Ir8a'}];
% Gen2Cons = [1:3];tic;
% [allMedBestFit,allPci_95,allPci_99] = getBestParameters(f_orco,stateKinAll,...
%     stateKinBaseline,Gen2Cons,lab,labCov,RegionLabel,10,false);
% save([meta.summationModelFold 'OrcoIr8aCombination.mat'])
% fprintf('Orco+Ir8a done in %d seconds.\n',round(toc));
% 
% %% Or42b Or92a
% lab = [{'Or42b','Or92a','Or42bOr92a','Or42b','Or92a'};
%  {'Before','Before','Before','',''}];
% 
% labCov = [{'Or42b'};
%     {'Or92a'}];
% Gen2Cons = [5,6,9];tic;
% [allMedBestFit,allPci_95,allPci_99] = getBestParameters(f_orco,stateKinAll,...
%     stateKinBaseline,Gen2Cons,lab,labCov,RegionLabel,10,false);
% save([meta.summationModelFold 'Or42bOr92aCombination.mat'])
% fprintf('Or42b+Or92a done in %d seconds.\n',round(toc));
% 
% %% Ir64a Ir75a
% lab = [{'Ir64a','Ir75a','Ir64aIr75a','Ir64a','Ir75a'};
%  {'Before','Before','Before','',''}];
% 
% labCov = [{'Ir64a'};
%     {'Ir75a'}];
% Gen2Cons = [7,8,10];tic;
% [allMedBestFit,allPci_95,allPci_99] = getBestParameters(f_orco,stateKinAll,...
%     stateKinBaseline,Gen2Cons,lab,labCov,RegionLabel,10,false);
% save([meta.summationModelFold 'Ir64aIr75aCombination.mat'])
% fprintf('Ir64a+Ir75a done in %d seconds.\n',round(toc));
% 
% %% Or42b Ir64a
% lab = [{'Or42b','Ir64a','Or42bIr64a','Or42b','Ir64a'};
%  {'Before','Before','Before','',''}];
% 
% labCov = [{'Or42b'};
%     {'Ir64a'}];
% Gen2Cons = [5,7,11];tic;
% [allMedBestFit,allPci_95,allPci_99] = getBestParameters(f_orco,stateKinAll,...
%     stateKinBaseline,Gen2Cons,lab,labCov,RegionLabel,10,false);
% save([meta.summationModelFold 'Or42bIr64aCombination.mat'])
% fprintf('Or42b+Ir64a done in %d seconds.\n',round(toc));
% 
% %% Ir64aIr75a Or42b
% lab = [{'Or42b','Ir64aIr75a','Ir64aIr75aOr42b','Or42b','Ir64aIr75a'};
%  {'Before','Before','Before','',''}];
% 
% labCov = [{'Or42b'};
%     {'Ir64aIr75a'}];
% Gen2Cons = [5,10,12];tic;
% [allMedBestFit,allPci_95,allPci_99] = getBestParameters(f_orco,stateKinAll,...
%     stateKinBaseline,Gen2Cons,lab,labCov,RegionLabel,10,false);
% save([meta.summationModelFold 'Ir64aIr75aOr42bCombination.mat'])
% fprintf('Or42b+Ir64aIr75a done in %d seconds.\n',round(toc));
% 
% %% Ir64aIr75a Ir75a
% lab = [{'Ir75a','Ir64aOr42b','Ir64aIr75aOr42b','Ir75a','Ir64aOr42b'};
%  {'Before','Before','Before','',''}];
% 
% labCov = [{'Ir75a'};
%     {'Ir64aOr42b'}];
% Gen2Cons = [8,11,12];tic;
% [allMedBestFit,allPci_95,allPci_99] = getBestParameters(f_orco,stateKinAll,...
%     stateKinBaseline,Gen2Cons,lab,labCov,RegionLabel,20,false);
% save([meta.summationModelFold 'Ir64aIr75aOr42bCombination2.mat'])
% fprintf('Ir64aOr42b+Ir75a done in %d seconds.\n',round(toc));
% 
% %% Or42a Or42b Or92a
% lab = [{'Or42a','Or42bOr92a','Or42aOr42bOr92a','Or42a','Or42bOr92a'};
%  {'Before','Before','Before','',''}];
% 
% labCov = [{'Or42a'};
%     {'Or42bOr92a'}];
% Gen2Cons = [4,9,13];tic;
% [allMedBestFit,allPci_95,allPci_99] = getBestParameters(f_orco,stateKinAll,...
%     stateKinBaseline,Gen2Cons,lab,labCov,RegionLabel,10,false);
% save([meta.summationModelFold 'Or42aOr42bOr92aCombination.mat'])
% fprintf('Or42a+Or42bOr92a done in %d seconds.\n',round(toc));

%% Plotting functions
allMat = {'OrcoIr8aCombination.mat','Or42bOr92aCombination.mat',...
    'Ir64aIr75aCombination.mat','Or42bIr64aCombination.mat',...
    'Ir64aIr75aOr42bCombination.mat','Ir64aIr75aOr42bCombination2.mat',...
    'Or42aOr42bOr92aCombination.mat'};
for i = 1:numel(allMat)
    allMat{i} = [meta.summationModelFold allMat{i}];
end

%plot posterior
close all
plotPosterior(allMat);

% plot synergy
close all
fNum = 1;
fNum = plotSynergy(allMat,'alpha',true,fNum);
fNum = plotSynergy(allMat,'full',true,fNum);
fNum = plotSynergy(allMat,'alpha',false,fNum);
fNum = plotSynergy(allMat,'full',false,fNum);
fNum = plotCovariance(allMat,fNum);

nORN = [36,5,41,1,1,1,1,1,2,2,2,3,3];
states2Cons = [1,2,3,4,7];
load(allMat{end},'stateKinBaseline','stateKinAll','RegionLabel','f_orco');
fNum = plotSynergisticByNumberOfORN(f_orco,stateKinBaseline,stateKinAll,...
    RegionLabel,nORN,states2Cons,fNum);

% states2Cons = [1,2,3,4];
% newOrder = [2,3,1,4,5];
% plotLooseRigid(allMat,states2Cons,newOrder,fNum);
% printFigures('Figures/ORN optogenetics_Synergy_Cov')
% 
% close all
% plotSyn = true;
% type = 'full';
% plotSynergyPosteriorByROI(allMat,type,plotSyn);
% printFigures('Figures/ORN optogenetics_SynergyPosteriorByROI_full')
% 
% close all
% plotSyn = true;
% type = 'alpha';
% plotSynergyPosteriorByROI(allMat,type,plotSyn);
% printFigures('Figures/ORN optogenetics_SynergyPosteriorByROI_alpha')
% 
% 
% close all
% type = 'absolute';%relative
% threshold = [0.5,0.5,1,0.5];
% plotSignificanceChangeByROI(allMat,t1*2.5,t2*2.5,t3,4.17+5,type,threshold);
% printFigures('Figures/ORN optogenetics_SynergyPosteriorByROI_realSpace_relative3')

close all
type = 'relative';%absolute
%threshold = [0.15,0.15,0.15,0.15];
threshold = [0.1,0.1,0.1,0.1];
plotSignificanceChangeByROI(allMat,t1*2.5,t2*2.5,t3,4.17+5,type,threshold);
printFigures('Figures/ORN optogenetics_SynergyPosteriorByROI_realSpace_relative3')

end

function [] = plotSynergyPosteriorByROI(allMat,type,plotSyn)

SynergyIdx = 11;
MuIdx = [4:5];

%plot covariance
nGen2Cons = 3;

pBaseline = zeros(nGen2Cons);
% for i = 1:nGen2Cons
%     pBaseline(i,i) = 1;
% end
pGen(1,:)  = [1,0];
pGen(2,:)  = [0,1];
pGen(3,:)  = [1,1];
pGen = [pBaseline,pGen];
pCovGen(:,1) = [0;0;1];

x = [2,8,4,6,5];

dFType = {'-df','0df','+df'};
fType = {'highF','lowF','Inh'};

states2Cons = [1,2,3,4];
for i = 1:numel(allMat)
    figure;set(gcf,'Position',[2 42 838 924])
    load(allMat{i},'allMedBestFit','f_orco','RegionLabel','labCov')
    if strcmpi(type,'alpha')
        colLab = {['mu ' labCov{1}],['mu ' labCov{2}],'muAB','Synergy (alpha)'};
    else
        colLab = {['mu ' labCov{1}],['mu ' labCov{2}],'muAB','Synergy (full)'};
    end
    covar = zeros(numel(states2Cons),5);
    
    varAB = cell(numel(states2Cons),4);
    synergy = cell(numel(states2Cons),4);
    k = 1;
    for state = states2Cons
        varAB(state,:) = {zeros(3,3)};
        synergy(state,:) = {zeros(3,3)};
        m = f_orco.model.params{state};
        stateType{1,state} = m.state.state;
        stateType{2,state} = m.state.kin;
        for ROI = 1:5
            tmp = allMedBestFit{state,ROI}(1:5)*pGen'+allMedBestFit{state,ROI}(11)*pCovGen';
            for j = 1:3
                varAB{state,j}(x(ROI)) = tmp(j);
            end
            if strcmpi(type,'alpha')
                synergy{state,4}(x(ROI)) = allMedBestFit{state,ROI}(SynergyIdx)./prod(allMedBestFit{state,ROI}(MuIdx));
            else
                synergy{state,4}(x(ROI)) = allMedBestFit{state,ROI}(SynergyIdx);
            end
            if plotSyn
                prodMat = sign(prod(allMedBestFit{state,ROI}(MuIdx)));
                synMat = sign((allMedBestFit{state,ROI}(SynergyIdx)));
                singleMat = sign(allMedBestFit{state,ROI}(MuIdx(1)));
                
                if prodMat==1 && synMat==singleMat
                    varAB{state,4}(x(ROI)) = 1;
                elseif prodMat==1 && synMat~=singleMat
                    varAB{state,4}(x(ROI)) = -1;
                else
                    varAB{state,4}(x(ROI)) = 0.01;
                end
            end
        end
        
        for j = 1:4
            subplot(numel(states2Cons),4,k)
            imagesc(varAB{state,j},[-1 1])
            
            set(gca,'ytick',[1:1:8],'yticklabel',fType);
            set(gca,'xtick',[1:1:5],'xticklabel',dFType);
            cMap = jet(201);
            cMap(101,:)=1;
            colormap(cMap);
            
            if j <4
                t2 = compose('%g',round(varAB{state,j},2));
                neg = varAB{state,j}<0;
                pos = varAB{state,j}>0;
            else
                t2 = compose('%g',round(synergy{state,j},2));
                neg = synergy{state,j}<0;
                pos = synergy{state,j}>0;
            end
            
            [xx,yy] = meshgrid([1:3],[1:3]);
            %text(xx(pos), yy(pos), t2(pos), 'HorizontalAlignment', 'Center','Color','k')%
            %text(xx(pos), yy(pos), t2(pos), 'HorizontalAlignment', 'Center','Color','k')%
            text(xx(:), yy(:), t2(:), 'HorizontalAlignment', 'Center','Color',0.5.*[1 1 1])
            title({[stateType{1,state} ' ' stateType{2,state}];colLab{j}})
            k = k+1;
        end
        cb=colorbar;cb.Position = cb.Position + [0.75e-1, 0, 0, 0];
        cb.YTick = [-1 0 1];
        cb.YTickLabel = {'Ant', 'N/A', 'Syn'};
    end
    sgtitle([labCov{1,1} ' + ' labCov{2,1}])
end

end

function [fNum] = plotLooseRigid(allMat,states2Cons,newOrder,fNum)

%plot covariance
nGen2Cons = 3;

pBaseline = zeros(nGen2Cons);
for i = 1:nGen2Cons
    pBaseline(i,i) = 1;
end

pGen(1,:)  = [1,0];
pGen(2,:)  = [0,1];
pGen(3,:)  = [1,1];
pGen = [pBaseline,pGen];

pCovGen(:,1) = [0;0;1];

for i = 1:numel(allMat)
    
    if mod(i,6)==1
        figure(fNum);set(gcf,'Position',[2 42 838 924]);
        fNum = fNum+1;
        k = 1;
    end
    
    %close all
    load(allMat{i},'allMedBestFit','f_orco','RegionLabel','labCov')
    covar = zeros(numel(states2Cons),5);
    for state = states2Cons
        m = f_orco.model.params{state};
        stateType{1,state} = m.state.state;
        stateType{2,state} = m.state.kin;
        for ROI = 1:5
            varAB = allMedBestFit{state,ROI}(6:10)*pGen'+2*allMedBestFit{state,ROI}(12)*pCovGen';
            if varAB(3)<min(varAB(1:2))
                covar(state,ROI) = -1;%rigid
            elseif varAB(3)>max(varAB(1:2))
                covar(state,ROI) = 1;%loose
            else
                covar(state,ROI) = 0;
            end
        end
    end
    
    subplot(3,2,k);
    imagesc(covar(:,newOrder),[-3 3]./2);
    tickLabels = strtrim(sprintf('%s\\newline%s\n', stateType{:,states2Cons}));
    set(gca,'ytick',[1:1:8],'yticklabel',tickLabels);
    set(gca,'xtick',[1:1:5],'xticklabel',RegionLabel(newOrder));
    title(strtrim(sprintf('%s + %s\n', labCov{:})));
    colormap(jet);
    if mod(k,2)==0 || i==numel(allMat)
        cb=colorbar;cb.Position = cb.Position + [0.75e-1, 0, 0, 0];
        cb.YTick = [-1 0 1];
        cb.YTickLabel = {'Rigid', 'N/A', 'Loose'};
    end
    
    k = k+1;
end

end

function [fNum] = plotSynergisticByNumberOfORN(f_orco,stateKinBaseline,...
    stateKinAll,RegionLabel,nORN,states2Cons,fNum)

nGen = size(stateKinAll,3);

figure(fNum);set(gcf,'Position',[2 42 838 924])
for stateNdx = 1:numel(states2Cons)
    state = states2Cons(stateNdx);
    m = f_orco.model.params{state};
    currBaselineMu = nanmean(log(stateKinBaseline{state,stateNdx}'));
    for ROI = 1:5
        allParams = nan(2,nGen);
        for i = 1:nGen
            tmpDat = log(stateKinAll{state,ROI,i}.allDat);
            mean2 = mean(tmpDat)-currBaselineMu;
            allParams(1,i) = mean2;
        end
        
        x0 = log10(nORN(1,:));
        y0 = allParams(1,:);
        p = polyfit(x0(~isnan(y0)),y0(~isnan(y0)),3);
        x1 = log10([1:41]);
        f1(1,:) = polyval(p,x1);
        
        subplot(numel(states2Cons),5,(stateNdx-1).*5+ROI);hold on;
        scatter(log10(nORN(1,:)),allParams(1,:),30,'r','LineWidth',1.5);
        plot(x1,f1(1,:),'-r','Linewidth',2)
        
        currMean = nanmean(allParams(1,nORN==1)).*[1:100];
        currUB = max(allParams(1,nORN==1)).*[1:100];
        currLB = min(allParams(1,nORN==1)).*[1:100];
        shadedErrorBar(log10(1:100),currMean,...
            [currUB-currMean;currMean-currLB])
        xlabel('log10(ORN)');
        if stateNdx == 1
            title(RegionLabel{ROI})
        end
        ylabel('mean');
        if ROI == 1
            ylabel({m.state.state ; m.state.kin})
        end
    end
end
sgtitle('#ORN to Mean, empirical distribution')
fNum = fNum+1;
end

function [fNum] = plotSynergy(allMat,type,plotSyn,fNum)
%plotSynergy
SynergyIdx = 11;
MuIdx = [4:5];
newOrder = [2,3,1,4,5];
tit = ['Synergy (' type ') between pairs'];

for i = 1:numel(allMat)
    %close all
    load(allMat{i})
    
    if mod(i,6)==1
        figure(fNum);set(gcf,'Position',[2 42 838 924]);
        suptitle(tit)
        fNum = fNum+1;
        k = 1;
    end
    
    synergy = nan(4,5);
    for state = 1:4
        m = f_orco.model.params{state};
        stateType{1,state} = m.state.state;
        stateType{2,state} = m.state.kin;
        for ROI = 1:5
            if strcmpi(type,'alpha')
                synergy(state,ROI) = allMedBestFit{state,ROI}(SynergyIdx)./prod(allMedBestFit{state,ROI}(MuIdx));
            else
                synergy(state,ROI) = allMedBestFit{state,ROI}(SynergyIdx);
            end
            if plotSyn
                prodMat(state,ROI) = sign(prod(allMedBestFit{state,ROI}(MuIdx)));
                synMat(state,ROI) = sign((allMedBestFit{state,ROI}(SynergyIdx)));
                singleMat(state,ROI) = sign(allMedBestFit{state,ROI}(MuIdx(1)));
            end
        end
    end
    
    subplot(3,2,k)
    if plotSyn==true
        tmpSynergy = nan(size(synergy));
        tmpSynergy(prodMat==1 & synMat==singleMat)=1;
        tmpSynergy(prodMat==1 & synMat~=singleMat)=-1;
        tmpSynergy(isnan(tmpSynergy)) = 0;
        imagesc(tmpSynergy(:,newOrder));
    else
        imagesc(synergy(:,newOrder));
    end
    tickLabels = strtrim(sprintf('%s\\newline%s\n', stateType{:}));
    set(gca,'ytick',[1:1:8],'yticklabel',tickLabels);
    set(gca,'xtick',[1:1:5],'xticklabel',RegionLabel(newOrder));
    colormap(jet);
    if mod(k,2)==0 || i==numel(allMat)
        cb=colorbar;cb.Position = cb.Position + [0.75e-1, 0, 0, 0];
        if plotSyn==true
            cb.YTick = [-1 0 1];
            cb.YTickLabel = {'Ant', 'N/A', 'Syn'};
        end
    end
    
    t2 = compose('%g',round(synergy(1:4,newOrder),2));
    neg = synergy(1:4,newOrder)<0;
    pos = synergy(1:4,newOrder)>0;
    
    [xx,yy] = meshgrid([1:5],[1:4]);
%     text(xx(pos), yy(pos), t2(pos), 'HorizontalAlignment', 'Center','Color','k')%
%     text(xx(neg), yy(neg), t2(neg), 'HorizontalAlignment', 'Center','Color','w')
    %text(xx(:), yy(:), t2, 'HorizontalAlignment', 'Center')%,'Color','w'
    text(xx(:), yy(:), t2(:), 'HorizontalAlignment', 'Center','Color',0.5.*[1 1 1])
    title([labCov{1,1} ' + ' labCov{2,1}])
    
    k = k+1;
end
end

function [fNum] = plotCovariance(allMat,fNum)
% plot covariance
CovarIdx = 12;
newOrder = [2,3,1,4,5];
tit = 'Covariance between pairs';


figure;set(gcf,'Position',[2 42 838 924])
for i = 1:numel(allMat)
    %close all
    load(allMat{i})
    
    if mod(i,6)==1
        figure(fNum);set(gcf,'Position',[2 42 838 924]);
        suptitle(tit)
        fNum = fNum+1;
        k = 1;
    end
    
    covar = nan(4,5);
    for state = 1:4
        m = f_orco.model.params{state};
        stateType{1,state} = m.state.state;
        stateType{2,state} = m.state.kin;
        for ROI = 1:5
            covar(state,ROI) = allMedBestFit{state,ROI}(CovarIdx)./2;
        end
    end
    subplot(3,2,k)
    imagesc(covar(:,newOrder));
    tickLabels = strtrim(sprintf('%s\\newline%s\n', stateType{:}));
    set(gca,'ytick',[1:1:8],'yticklabel',tickLabels);
    set(gca,'xtick',[1:1:5],'xticklabel',RegionLabel(newOrder));
    colormap(jet);
    if mod(k,2)==0 || i==numel(allMat)
        cb=colorbar;cb.Position = cb.Position + [0.75e-1, 0, 0, 0];
    end
    
    t2 = compose('%g',round(covar(1:4,newOrder),2));
    neg = covar(1:4,newOrder)<0;
    pos = covar(1:4,newOrder)>0;
    
    [xx,yy] = meshgrid([1:5],[1:4]);
%     text(xx(pos), yy(pos), t2(pos), 'HorizontalAlignment', 'Center','Color','k')%
%     text(xx(neg), yy(neg), t2(neg), 'HorizontalAlignment', 'Center','Color','w')
    text(xx(:), yy(:), t2(:), 'HorizontalAlignment', 'Center','Color',0.5.*[1 1 1])
    %text(xx(:), yy(:), t2, 'HorizontalAlignment', 'Center')%,'Color','w'
    title([labCov{1,1} ' + ' labCov{2,1}])
    
    k = k+1;
end
end

function [] = plotPosterior(allMat)
for i = 1:numel(allMat)
    close all
    load(allMat{i})
    %%
    plotDistributionFit(f_orco,allMedBestFit,stateKinAll,RegionLabel,genAll,Gen2Cons);
    %%
    plotPosteriorMu(f_orco,allMedBestFit,allPci_95,allPci_99,RegionLabel,lab);
    %%
    plotPosteriorVar(f_orco,allMedBestFit,allPci_95,allPci_99,RegionLabel,lab);
    %%
    plotPosteriorDistribution(f_orco,allMedBestFit,RegionLabel,lab);
    
    C = strsplit(allMat{i},{'.','\'});
    fName = [C{end-1} '_Posterior'];
    %fName = [allMat{i}(1:end-4) '_Posterior'];
    printFigures(['Figures/' fName]);
end
end

function [] = plotPosteriorDistribution(f_orco,allMedBestFit,RegionLabel,lab)
k = 1;
nParam = 5;
figure;set(gcf,'Position',[2 42 838 924])
for state = 1:8
    m = f_orco.model.params{state};
    for ROI = 1:5
        medBestFit = allMedBestFit{state,ROI};
        if ~isempty(medBestFit)
            
            medBestFitMu = medBestFit(1:nParam);
            medBestFitSigma = sqrt(abs(medBestFit(nParam+1:nParam*2)));
            sgn = sign(medBestFit(nParam+1:nParam*2));
            
            subplot(8,5,k);hold on;
            scatter(medBestFitMu,1:5,'ok')
            
            for i = 1:nParam
                if sgn(i)==1
                    plot(medBestFitMu(i)+2.*medBestFitSigma(i).*[-1 1],[i i],'k');
                else
                    plot(medBestFitMu(i)+2.*medBestFitSigma(i).*[-1 1],[i i],'r');
                end
            end
            xlim([-5 5])
            if ROI == 1
                tickLabels = strtrim(sprintf('%s %s\n', lab{:}));
                set(gca,'ytick',[0.5:1:5.5],'yticklabel',tickLabels);
                xlabel([m.state.state ' ' m.state.kin]);
            end
            if state == 1
                title(RegionLabel{ROI})
            end
        end
        k = k+1;
    end
end
sgtitle(['ORN contribution, mean+/-2|STD|, red=(-) variance'])
end

function [] = plotPosteriorVar(f_orco,allMedBestFit,allPci_95,allPci_99,RegionLabel,lab)
k = 1;
nParam = 5;
figure;set(gcf,'Position',[2 42 838 924])
for state = 1:8
    m = f_orco.model.params{state};
    for ROI = 1:5
        medBestFit = allMedBestFit{state,ROI};
        if ~isempty(medBestFit)
            pci_95 = allPci_95{state,ROI};
            pci_99 = allPci_99{state,ROI};
            
            medBestFitVar = medBestFit(nParam+1:2*nParam);
            pciVar_95 = pci_95(:,nParam+1:2*nParam);
            pciVar_99 = pci_99(:,nParam+1:2*nParam);
            
            subplot(8,5,k);hold on;
            for i = 1:nParam
                plot(pciVar_95(:,i),[i i],'k','LineWidth',2.5);
                plot(pciVar_99(:,i),[i i],'k');
            end
            scatter(medBestFitVar,1:5,'ok')
            
            if ROI == 1
                tickLabels = strtrim(sprintf('%s %s\n', lab{:}));
                set(gca,'ytick',[0.5:1:5.5],'yticklabel',tickLabels);
                xlabel([m.state.state ' ' m.state.kin]);
            end
            if state == 1
                title(RegionLabel{ROI})
            end
        end
        k = k+1;
    end
end
sgtitle(['Variance MLE confidence, 95% and 99%'])
end

function [] = plotPosteriorMu(f_orco,allMedBestFit,allPci_95,allPci_99,RegionLabel,lab)
k = 1;
nParam = 5;
figure;set(gcf,'Position',[2 42 838 924])
for state = 1:8
    m = f_orco.model.params{state};
    for ROI = 1:5
        medBestFit = allMedBestFit{state,ROI};
        if ~isempty(medBestFit)
            pci_95 = allPci_95{state,ROI};
            pci_99 = allPci_99{state,ROI};
            
            medBestFitMu = medBestFit(1:nParam);
            pciMu_95 = pci_95(:,1:nParam);
            pciMu_99 = pci_99(:,1:nParam);
            
            subplot(8,5,k);hold on;
            for i = 1:nParam
                plot(pciMu_95(:,i),[i i],'k','LineWidth',2.5);
                plot(pciMu_99(:,i),[i i],'k');
            end
            scatter(medBestFitMu,1:5,'ok')
            
            if ROI == 1
                tickLabels = strtrim(sprintf('%s %s\n', lab{:}));
                set(gca,'ytick',[0.5:1:5.5],'yticklabel',tickLabels);
                xlabel([m.state.state ' ' m.state.kin]);
            end
            if state == 1
                title(RegionLabel{ROI})
            end
        end
        k = k+1;
    end
end
sgtitle(['Mu MLE confidence, 95% and 99%'])
end

function [] = plotDistributionFit(f_orco,allMedBestFit,stateKinAll,RegionLabel,genAll,Gen2Cons)

nGen2Cons = 3;
pBaseline = zeros(nGen2Cons);
for i = 1:nGen2Cons
    pBaseline(i,i) = 1;
end
pGen(1,:)  = [1,0];
pGen(2,:)  = [0,1];
pGen(3,:)  = [1,1];
pGen = [pBaseline,pGen];

pCovGen(1,:) = [0];
pCovGen(2,:) = [0];
pCovGen(3,:) = [1];
nParam = 5;
dx = [{-6:0.5:6};{-6:0.5:6};{0:0.5:6};{-6:0.5:6};{-6:0.5:6};{-6:0.5:6};{-10:0.5:10};{-10:0.5:10}];


for i = 1:3
    figure;set(gcf,'Position',[2 42 838 924])
    for state = 1:8
        m = f_orco.model.params{state};
        for ROI = 1:5
            if numel(allMedBestFit{state,ROI})>0
                currMu = allMedBestFit{state,ROI}(1:nParam)*pGen';
                currSig = allMedBestFit{state,ROI}(nParam+1:2*nParam)*pGen';
                currMuCov = allMedBestFit{state,ROI}(2*nParam+1:2*nParam+1)*pCovGen';
                currSigCov = allMedBestFit{state,ROI}(2*nParam+2:2*nParam+2)*pCovGen';
                tmpDat = log(stateKinAll{state,ROI,Gen2Cons(i)}.allDat);
                
                currMu = currMu(i)+currMuCov(i);
                currSig = sqrt(currSig(i)+2.*currSigCov(i));
                
                dxOpt = getOptBinSize(tmpDat,'Sturges');
                x_values = dx{state};
                y1 = normpdf(x_values,currMu,currSig);
                
                if ~isempty(tmpDat) && dxOpt>0
                    subplot(8,5,(state-1).*5+ROI);
                    histogram(tmpDat,x_values,'Normalization','pdf');hold on;
                    plot(x_values,y1,'k','Linewidth',2)
                    ylabel('pdf')
                    xlabel([m.state.state ' ' m.state.kin]);
                    if state == 1
                        title(RegionLabel{ROI})
                    end
                    set(gca,'xtick',[-6:2:6]);
                end
            end
        end
    end
    sgtitle([genAll{Gen2Cons(i)} ' bar=log(emp), line=fit'])
end
end

function [allMedBestFit,allPci_95,allPci_99] = getBestParameters...
    (f_orco,stateKinAll,stateKinBaseline,Gen2Cons,lab,labCov,RegionLabel,nIt,plotfig)
warning('off','stats:mle:IterLimit')
warning('off','stats:mle:EvalLimit')
nGen2Cons = 3;
nParam = 2+nGen2Cons;
nParamCov = size(labCov,2);

pBaseline = zeros(nGen2Cons);
pCovBaseline = zeros(3,1);
for i = 1:nGen2Cons
    pBaseline(i,i) = 1;
end

pGen(1,:)  = [1,0];
pGen(2,:)  = [0,1];
pGen(3,:)  = [1,1];
pGen = [pBaseline,pGen];
pBaseline = [pBaseline,zeros(size(pGen,1),size(pGen,2)-nGen2Cons)];

pCovGen(1,:) = [0];
pCovGen(2,:) = [0];
pCovGen(3,:) = [1];

tic;%nIt = 5;
for state = 1:8
    m = f_orco.model.params{state};
    
    if plotfig
        figure;set(gcf,'Position',[2 42 1069 924])
    end
    for ROI = 1:5
        
        y = [];pAll = [];pCovAll = [];
        yBaseline = [];pAllBaseline = [];pCovAllBaseline = [];
        for i = 1:numel(Gen2Cons)
            nPt = numel(stateKinAll{state,ROI,Gen2Cons(i)}.allDat);
            if nPt>15
                yBaseline = [yBaseline; stateKinBaseline{state,Gen2Cons(i)}];
                pAllBaseline = [pAllBaseline; repmat(pBaseline(i,:),numel(stateKinBaseline{state,Gen2Cons(i)}),1)];
                pCovAllBaseline = [pCovAllBaseline; repmat(pCovBaseline(i,:),numel(stateKinBaseline{state,Gen2Cons(i)}),1)];
                
                y = [y; stateKinAll{state,ROI,Gen2Cons(i)}.allDat];
                pAll = [pAll; repmat(pGen(i,:),numel(stateKinAll{state,ROI,Gen2Cons(i)}.allDat),1)];
                pCovAll = [pCovAll; repmat(pCovGen(i,:),numel(stateKinAll{state,ROI,Gen2Cons(i)}.allDat),1)];
            else
                a = 1;
            end
        end
        y = log([y;yBaseline]);
        pAll = [pAll;pAllBaseline];
        pCovAll = [pCovAll;pCovAllBaseline];
        
        if ~isempty(y)
            
            muF = @(p_mu,p_synergy) (p_mu*pAll'+p_synergy*pCovAll')';
            varF = @(p_sig,p_covar) (p_sig*pAll'+2.*p_covar*pCovAll')';
            pdf = @(y,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22) ...
                (1./(sqrt(varF([p7,p8,p9,p10,p11,p12],[p18,p19,p20,p21,p22])).*sqrt(2*pi)) .* ...
                exp(-((y-muF([p1,p2,p3,p4,p5,p6],[p13,p14,p15,p16,p17])).^2 ./ (2*(varF([p7,p8,p9,p10,p11,p12],[p18,p19,p20,p21,p22]))))))+eps;%+eps;
            nLL = @(y,p) -sum(log(pdf(y,p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10),p(11),p(12),p(13),p(14),p(15),p(16),p(17),p(18),p(19),p(20),p(21),p(22))+eps));
            %lb = [-5.*ones(1,nParam),0.*ones(1,nParam),-5.*ones(1,nParamCov.*2)];
            %ub = [6.*ones(1,nParam.*2+nParamCov.*2)];
            
            bestFit = zeros(nIt,nParam.*2+nParamCov.*2);
            pci_99= cell(nIt,1);pci_95= cell(nIt,1);
            for it = 1:nIt
                p_mu = rand(1,nParam).*4;
                p_var = rand(1,nParam).*1+1;
                p_synergy = 0.1.*rand(1,nParamCov).*1-0.5;
                p_covar = rand(1,nParamCov).*1-0.5;
                opt = statset('mlecustom');
                %opt.MaxIter = 1000;
                
                %tic;
                try
                    phat2(1,:) = mle(y,'pdf',@(data,varargin) CustumLogNormalPDF(data,pAll,pCovAll,5,1,varargin),...
                        'start',[p_mu p_var p_synergy p_covar],'Options',opt,'OptimFun', 'fminsearch');
                    nLL_all(1) = CustumLogNormalnLL(y,pAll,pCovAll,5,1,phat2(1,:));
                    %fprintf('iteration %d took %d seconds has a LL of %d.\n',i,round(toc),LL(i));
                    for j = 2:200
                        phat2(j,:) = mle(y,'pdf',@(data,varargin) CustumLogNormalPDF(data,pAll,pCovAll,5,1,varargin),...
                            'start',phat2(j-1,:),'Options',opt,'OptimFun', 'fminsearch');
                        nLL_all(j) = CustumLogNormalnLL(y,pAll,pCovAll,5,1,phat2(j,:));
                        if abs(nLL_all(j)-nLL_all(j-1))<1
                            break
                        end
                    end
                    if nLL_all(j)<=nLL_all(j-1)
                        [phat2(j,:), pci_95{it}] = mle(y,'pdf',@(data,varargin) CustumLogNormalPDF(data,pAll,pCovAll,5,1,varargin),...
                            'start',phat2(j-1,:),'Options',opt,'OptimFun', 'fmincon');
                        [phat2(j,:), pci_99{it}] = mle(y,'pdf',@(data,varargin) CustumLogNormalPDF(data,pAll,pCovAll,5,1,varargin),...
                            'start',phat2(j-1,:),'Options',opt,'OptimFun', 'fmincon','Alpha',0.01);
                        bestFit(it,:) = phat2(j,:);
                    else
                        [phat2(j-1,:), pci_95{it}] = mle(y,'pdf',@(data,varargin) CustumLogNormalPDF(data,pAll,pCovAll,5,1,varargin),...
                            'start',phat2(j-2,:),'Options',opt,'OptimFun', 'fmincon');
                        [phat2(j-1,:), pci_99{it}] = mle(y,'pdf',@(data,varargin) CustumLogNormalPDF(data,pAll,pCovAll,5,1,varargin),...
                            'start',phat2(j-2,:),'Options',opt,'OptimFun', 'fmincon','Alpha',0.01);
                        bestFit(it,:) = phat2(j-1,:);
                    end
                catch
                    bestFit(it,:) = nan(1,nParam.*2+nParamCov.*2);
                end
%                 try
%                     fprintf('iteration %d took %d seconds has a LL of %d.\n',it,round(toc),nLL_all(end));
%                 catch
%                     a = 1;
%                 end
            end
            pci_99(isnan(bestFit(:,1))) = [];
            pci_95(isnan(bestFit(:,1))) = [];
            bestFit(isnan(bestFit(:,1)),:) = [];
            
            badNdx = cellfun(@(x) any(isnan(x),'all'),pci_95)==1 | cellfun(@(x) any(isnan(x),'all'),pci_99)==1;
            if ~all(badNdx)
                pci_99(badNdx) = [];
                pci_95(badNdx) = [];
                bestFit(badNdx,:) = [];
            end
            
            medBestFit = median(bestFit,1);
            [~,BestFitNdx] = min(sum(abs(bestFit-medBestFit),2));
            %input = num2cell(medBestFit);
            %fitPDF = pdf(y,input{:});
            
            if plotfig
                for i = 1:nParam
                    subplot(5,4,4.*(ROI-1)+1);hold on;
                    scatter(0.5.*linspace(0,1,size(bestFit,1))+i+eps,bestFit(:,i));
                    plot([i,i+0.5],[median(bestFit(:,i)), median(bestFit(:,i))],'k','linewidth',2)
                    xlim([0.5 nParam+1])
                    subplot(5,4,4.*(ROI-1)+2);hold on;
                    scatter(0.5.*linspace(0,1,size(bestFit,1))+i,bestFit(:,i+nParam));
                    plot([i,i+0.5],[median(bestFit(:,i+nParam)), median(bestFit(:,i+nParam))],'k','linewidth',2)
                    xlim([0.5 nParam+1])
                    
                    if i<nParamCov+1
                        subplot(5,4,4.*(ROI-1)+3);hold on;
                        scatter(0.5.*linspace(0,1,size(bestFit,1))+i,bestFit(:,i+2*nParam));
                        plot([i,i+0.5],[median(bestFit(:,i+2*nParam)), median(bestFit(:,i+2*nParam))],'k','linewidth',2)
                        xlim([0.5 nParamCov+1])
                        
                        subplot(5,4,4.*ROI);hold on;
                        scatter(0.5.*linspace(0,1,size(bestFit,1))+i,bestFit(:,i+2*nParam+nParamCov));
                        plot([i,i+0.5],[median(bestFit(:,i+2*nParam+nParamCov)), median(bestFit(:,i+2*nParam+nParamCov))],'k','linewidth',2)
                        xlim([0.5 nParamCov+1])
                    end
                end
                tickLabels = strtrim(sprintf('%s\\newline%s\n', lab{:}));
                subplot(5,4,4.*(ROI-1)+1);set(gca,'xtick',[1:nParam],'xticklabel',tickLabels)
                xtickangle(20)
                title([RegionLabel{ROI} ' mu']);ylabel('mu(log(x))')
                subplot(5,4,4.*(ROI-1)+2);set(gca,'xtick',[1:nParam],'xticklabel',tickLabels)
                xtickangle(20)
                title([RegionLabel{ROI} ' var']);ylabel('var(log(x))')
                
                tickLabels = strtrim(sprintf('%s\\newline%s\n', labCov{:}));
                subplot(5,4,4.*(ROI-1)+3);set(gca,'xtick',[1:nParamCov],'xticklabel',tickLabels)
                xtickangle(20)
                title([RegionLabel{ROI} ' syn']);ylabel('syn(log(x))')
                subplot(5,4,4.*ROI);set(gca,'xtick',[1:nParamCov],'xticklabel',tickLabels)
                xtickangle(20)
                title([RegionLabel{ROI} ' covar']);ylabel('covar(log(x))')
                suptitle([m.state.state ' ' m.state.kin]);
            end
            try
                allMedBestFit{state,ROI} = medBestFit;
                allPci_95{state,ROI} = pci_95{BestFitNdx};
                allPci_99{state,ROI} = pci_99{BestFitNdx};
            catch
                allMedBestFit{state,ROI} = [];
                allPci_95{state,ROI} = [];
                allPci_99{state,ROI} = [];
            end
        end
    end
    %toc;
    fprintf('state %d /8 done in %d seconds.\n',state,round(toc));
end
end

function [stateKinAll,stateKinBaseline,stateKinBaselineMu,stateKinBaselineVar,...
    StateLabel,f_orco] = getKinematicsByRegions(genAll,meta,t1,t2,t3,plotFig)
for i = 1 :numel(genAll)%2%
    gen = genAll{i};
    
    load(['DataModel/' gen '_' meta.d meta.ext '.mat'],'f_orco');
    baseline = f_orco.spk(1);
    
    XX = f_orco.model.TurnBias.XX;
    YY = f_orco.model.TurnBias.YY;
    for state = 1:8
        k = 1;
        m = f_orco.model.params{state};
        baseLineNdx = abs(m.alldSpk)<0.1 & (m.allSpk-baseline)<0.1;
        
        if (i == 1) && (state == 2) && plotFig
            figure;scatter(m.alldSpk(~baseLineNdx),m.allSpk(~baseLineNdx));hold on;
            plot([t1 t1],[0 50],'k','Linewidth',1);
            plot([t2 t2],[0 50],'k','Linewidth',1);
            plot([t1 t2],[baseline baseline],'k','Linewidth',1);
            plot([t1 t2],[t3 t3],'k','Linewidth',1);
            xlim([-200 200]);ylabel([0 45]);
            xlabel('df (spks/s^2)');ylabel('f (spks/s)')
            text(-150,40,'NegDF','FontWeight','Bold','FontSize',12,'HorizontalAlignment','Center');
            text(150,40,'PosDF','FontWeight','Bold','FontSize',12,'HorizontalAlignment','Center');
            text(0,40,'HighF','FontWeight','Bold','FontSize',12,'HorizontalAlignment','Center');
            text(0,2.5,'Inh','FontWeight','Bold','FontSize',12,'HorizontalAlignment','Center');
            text(0,10,'LowF','FontWeight','Bold','FontSize',12,'HorizontalAlignment','Center');
        end
        nonBaselineDat = m.allDat(~baseLineNdx)+10*eps;
        nonBaselineDSpk = m.alldSpk(~baseLineNdx);
        nonBaselineSpk = m.allSpk(~baseLineNdx);
        BaselineDat = m.baselineDat+10*eps;
        
        ROI1 = nonBaselineDSpk<t1;
        ROI2 = nonBaselineDSpk>t2;
        ROI3 = (nonBaselineSpk>t3) & ~ROI1 & ~ROI2;
        ROI4 = (nonBaselineSpk<baseline) & ~ROI1 & ~ROI2;
        ROI5 = ~ROI1 & ~ROI2 & ~ROI3 & ~ROI4;
        
        stateKinAll{state,1,i}.alldSpk = nonBaselineDSpk(ROI1);
        stateKinAll{state,2,i}.alldSpk = nonBaselineDSpk(ROI2);
        stateKinAll{state,3,i}.alldSpk = nonBaselineDSpk(ROI3);
        stateKinAll{state,4,i}.alldSpk = nonBaselineDSpk(ROI4);
        stateKinAll{state,5,i}.alldSpk = nonBaselineDSpk(ROI5);
        stateKinAll{state,1,i}.allSpk = nonBaselineSpk(ROI1);
        stateKinAll{state,2,i}.allSpk = nonBaselineSpk(ROI2);
        stateKinAll{state,3,i}.allSpk = nonBaselineSpk(ROI3);
        stateKinAll{state,4,i}.allSpk = nonBaselineSpk(ROI4);
        stateKinAll{state,5,i}.allSpk = nonBaselineSpk(ROI5);
        stateKinAll{state,1,i}.allDat = nonBaselineDat(ROI1);
        stateKinAll{state,2,i}.allDat = nonBaselineDat(ROI2);
        stateKinAll{state,3,i}.allDat = nonBaselineDat(ROI3);
        stateKinAll{state,4,i}.allDat = nonBaselineDat(ROI4);
        stateKinAll{state,5,i}.allDat = nonBaselineDat(ROI5);
        
        stateKinBaseline{state,i} = BaselineDat;
        stateKinBaselineMu(state,i) = mean(log(BaselineDat));
        stateKinBaselineVar(state,i) = var(log(BaselineDat));
        
        StateLabel{state} = [m.state.state ' ' m.state.kin];
    end
end
end

function [] = printFigures(fName)

for f = 1:get(gcf,'Number')
    figure(f);
    print('-painters','-dpsc2',[fName '.ps'],'-append');
end
ps2pdf('psfile', [fName '.ps'], 'pdffile', [fName '.pdf'], 'gspapersize', 'letter',...
    'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
    'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
    'gslibpath','C:\Program Files\gs\gs9.50\lib');


end


