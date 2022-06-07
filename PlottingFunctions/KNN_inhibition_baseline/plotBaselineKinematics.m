function [fNum] = plotBaselineKinematics(f_orcoAll,border,fNum)

fNumInitial = fNum;
fNum = fNum+2;

nGen = numel(f_orcoAll);
stateDistEarly = cell(4,nGen);stateDistLate = cell(4,nGen);stateDistBefore = cell(4,nGen);
baselineEarlyRankSum = nan(7,nGen);
baselineLateRankSum = nan(7,nGen);
baselineAllRankSum = nan(7,nGen);
baselineEarlyRankSum_increase = nan(7,nGen);
baselineEarlyRankSum_decrease = nan(7,nGen);
baselineLateRankSum_increase = nan(7,nGen);
baselineLateRankSum_decrease = nan(7,nGen);
baselineAllRankSum_increase = nan(7,nGen);
baselineAllRankSum_decrease = nan(7,nGen);
allGenName = cell(3,nGen);
nMin = 30;
for gen = 1:nGen
    f_orco = f_orcoAll{gen};
    tmpGenName = erase(f_orco.id,{' Retinal',' Control'});
    
    tmp=split(tmpGenName,' ');
    allGenName(1:numel(tmp),gen) = tmp;
    
    baseline = getBaselineKinematics2(f_orco,border,[],[]);%baseline{gen}
    
    allTracks = baseline.during;
    allStates = allTracks.state;

    earlyNdx = 1:2;
    lateNdx = earlyNdx(end)+1:10;%size(allStates,2);

    nStates = numel(f_orco.states.key);
    
    for state = 1:nStates
        stateDistEarly{state,gen}.curv = [];
        stateDistEarly{state,gen}.spd = [];
        stateDistEarly{state,gen}.dur = [];

        stateDistLate{state,gen}.curv = [];
        stateDistLate{state,gen}.spd = [];
        stateDistLate{state,gen}.dur = [];

        stateDistBefore{state,gen}.curv = [];
        stateDistBefore{state,gen}.spd = [];
        stateDistBefore{state,gen}.dur = [];
    end
    
    [stateDistBefore(:,gen)] = getKinematics(f_orco,baseline.before.state,baseline.before,stateDistBefore(:,gen),[1:1000]);
%     % sanity check
%     if ~isempty(f_orco.model.params)
%         % may contain 1 extra track in some flies (<1% difference)
%         tmp2 = [nanmean(stateDistBefore{1,gen}.spd), nanmean(f_orco.model.params{1,1}.baselineDat');...
%         nanmean(stateDistBefore{2,gen}.spd), nanmean(f_orco.model.params{1,2}.baselineDat');...
%         nanmean(stateDistBefore{1,gen}.curv), nanmean(f_orco.model.params{1,3}.baselineDat');...
%         nanmean(stateDistBefore{2,gen}.curv), nanmean(f_orco.model.params{1,4}.baselineDat');...
%         nanmean(stateDistBefore{1,gen}.dur), nanmean(f_orco.model.params{1,5}.baselineDat');...
%         nanmean(stateDistBefore{2,gen}.dur), nanmean(f_orco.model.params{1,6}.baselineDat');...
%         nanmean(stateDistBefore{3,gen}.curv), nanmean(f_orco.model.params{1,7}.baselineDat');...
%         nanmean(stateDistBefore{3,gen}.dur), nanmean(f_orco.model.params{1,8}.baselineDat')];
%     end
    
    [stateDistEarly(:,gen)] = getKinematics(f_orco,allStates,allTracks,stateDistEarly(:,gen),earlyNdx);
    [stateDistLate(:,gen)] = getKinematics(f_orco,allStates,allTracks,stateDistLate(:,gen),lateNdx);
    
    k = 1;
    for i = 1:3
        fn = fieldnames(stateDistBefore{i,gen});
        if i==3
            fn = fn(3);
        end
        for j = 1:numel(fn)
            stateType{1,k} = f_orco.states.key{i};
            stateType{2,k} = fn{j};
            
            if numel(stateDistEarly{i,gen}.(fn{j}))>nMin
                [baselineEarlyRankSum(k,gen),~] = ranksum(stateDistEarly...
                    {i,gen}.(fn{j}),stateDistBefore{i,gen}.(fn{j}));
                [baselineEarlyRankSum_increase(k,gen),~] = ranksum(stateDistEarly...
                    {i,gen}.(fn{j}),stateDistBefore{i,gen}.(fn{j}),'tail','right');
                [baselineEarlyRankSum_decrease(k,gen),~] = ranksum(stateDistEarly...
                    {i,gen}.(fn{j}),stateDistBefore{i,gen}.(fn{j}),'tail','left');
            end
            if numel(stateDistLate{i,gen}.(fn{j}))>nMin
                [baselineLateRankSum(k,gen),~] = ranksum(stateDistLate...
                    {i,gen}.(fn{j}),stateDistBefore{i,gen}.(fn{j}));
                [baselineLateRankSum_increase(k,gen),~] = ranksum(stateDistLate...
                    {i,gen}.(fn{j}),stateDistBefore{i,gen}.(fn{j}),'tail','right');
                [baselineLateRankSum_decrease(k,gen),~] = ranksum(stateDistLate...
                    {i,gen}.(fn{j}),stateDistBefore{i,gen}.(fn{j}),'tail','left');
            end
            if numel([stateDistEarly{i,gen}.(fn{j}),stateDistLate{i,gen}.(fn{j})])>nMin
                tmpAllDist = [stateDistEarly{i,gen}.(fn{j}),stateDistLate{i,gen}.(fn{j})];
                [baselineAllRankSum(k,gen),~] = ranksum(tmpAllDist,...
                    stateDistBefore{i,gen}.(fn{j}));
                [baselineAllRankSum_increase(k,gen),~] = ranksum(tmpAllDist,...
                    stateDistBefore{i,gen}.(fn{j}),'tail','right');
                [baselineAllRankSum_decrease(k,gen),~] = ranksum(tmpAllDist,...
                    stateDistBefore{i,gen}.(fn{j}),'tail','left');
            end
            nAll{k,gen} = [numel(stateDistBefore{i,gen}.(fn{j})),...
                numel(stateDistEarly{i,gen}.(fn{j})),numel(stateDistLate{i,gen}.(fn{j}))];
            k=k+1;
        end
    end
    
    plotDistributions(f_orco,stateDistBefore(:,gen),stateDistEarly(:,gen),stateDistLate(:,gen),earlyNdx,lateNdx,nMin,fNum);
    suptitle([f_orco.id ' baseline Tracks'])
    fNum = fNum+1;
end

pThresh = 0.05;
tmpAll = {baselineEarlyRankSum,baselineLateRankSum,baselineAllRankSum};
%{baselineAllRankSum_increase,baselineAllRankSum_decrease}
for i = 1:3
    if i == 1
        a = baselineEarlyRankSum_increase;
        b = baselineAllRankSum_decrease;
    elseif i == 2
        a = baselineLateRankSum_increase;
        b = baselineLateRankSum_decrease;
    elseif i == 3
        a = baselineAllRankSum_increase;
        b = baselineAllRankSum_decrease;
    end
    tmpAll{i} = min(cat(3,a,b),[],3);
    tmp{i} = zeros(size(a));
    tmp{i}(a<pThresh) = 1;
    tmp{i}(b<pThresh) = -1;
    tmp{i}(isnan(a)) = -2;
end
[xx,yy] = meshgrid([1:nGen],[1:7]);
cMap = [0,0,0;0,0,1;1,1,1;1,0,0];
figure(fNumInitial);set(gcf,'Position',[2 42 838 924])
for i = 1:3
    subplot(3,1,i);
    imagesc(tmp{i},[-2 1])
    kinTickLabels = strtrim(sprintf('%s\\newline%s\n', stateType{:}));
    genTickLabels = strtrim(sprintf('%s\\newline%s\\newline%s\n', allGenName{:}));
    set(gca,'ytick',[1:1:7],'yticklabel',kinTickLabels);
    set(gca,'xtick',[1:1:nGen],'xticklabel',genTickLabels);
    colormap(cMap);
    cb=colorbar;cb.Position = cb.Position + [0.75e-1, 0, 0, 0];
    cb.YTick = [-2 -1 0 1];
    cb.YTickLabel = {'N/A', 'Lower', 'NS', 'Higher'};
    
    t2 = ones(size(tmpAll{i}));
    t2(tmp{i}==1 | tmp{i}==-1) = tmpAll{i}(tmp{i}==1 | tmp{i}==-1);
    t2 = compose('%g',round(t2,3));
    t2(strcmpi(t2,{'1'}))={''};
    text(xx(:), yy(:), t2(:), 'HorizontalAlignment', 'Center','Color','w')
end
subplot(3,1,1);title('Rank Sum (one sided) Before vs first 2 baseline states')
subplot(3,1,2);title('Rank Sum (one sided) Before vs later baseline states')
subplot(3,1,3);title('Rank Sum (one sided) Before vs first all baseline states')

pThresh = 0.05;
tmpAll = {baselineEarlyRankSum,baselineLateRankSum,baselineAllRankSum};
for i = 1:3
    tmp{i} = zeros(size(tmpAll{i}));
    tmp{i}(tmpAll{i}<pThresh) = 1;
    tmp{i}(isnan(tmpAll{i})) = -1;
end
[xx,yy] = meshgrid([1:nGen],[1:7]);

cMap = ones(1,3).*[1;1;0];
cMap(1,:) = [1,0,0];
figure(fNumInitial+1);set(gcf,'Position',[2 42 838 924])
for i = 1:3
    subplot(3,1,i);
    imagesc(tmp{i})
    kinTickLabels = strtrim(sprintf('%s\\newline%s\n', stateType{:}));
    genTickLabels = strtrim(sprintf('%s\\newline%s\\newline%s\n', allGenName{:}));
    set(gca,'ytick',[1:1:7],'yticklabel',kinTickLabels);
    set(gca,'xtick',[1:1:nGen],'xticklabel',genTickLabels);
    colormap(cMap);
    cb=colorbar;cb.Position = cb.Position + [0.75e-1, 0, 0, 0];
    cb.YTick = [-1 0 1];
    cb.YTickLabel = {'N/A', 'NS', 'Sig'};
    
    t2 = ones(size(tmpAll{i}));
    t2(tmp{i}==1) = tmpAll{i}(tmp{i}==1);
    t2 = compose('%g',round(t2,3));
    t2(strcmpi(t2,{'1'}))={''};
    text(xx(:), yy(:), t2(:), 'HorizontalAlignment', 'Center','Color','w')
end
subplot(3,1,1);title('Rank Sum Before vs first 2 baseline states')
subplot(3,1,2);title('Rank Sum Before vs later baseline states')
subplot(3,1,3);title('Rank Sum Before vs first all baseline states')

end


function [stateDist] = getKinematics(self,allStates,allTracks,stateDist,ndx)

nStates = numel(self.states.key);
for fly = 1:self.nFly
    currFlyStates = squeeze(allStates(fly,:,:));
    goodTracks = find(~all(isnan(currFlyStates),2));
    
    for track = 1:numel(goodTracks)
        currTrackStates = currFlyStates(goodTracks(track),:);
        currTrackStates(isnan(currTrackStates)) = [];
        
        for state = 1:nStates
            [startNdx,endNdx,type] = startEndSeq(currTrackStates==state);
            startNdx = startNdx(type);endNdx = endNdx(type);
            len = endNdx-startNdx+1;
            startNdx(len<3) = [];endNdx(len<3) = [];
            
            for i = 1:numel(startNdx)
                if strcmpi(self.states.key{state},'sharp turns')
                    tmpCurve = abs(sum(allTracks.curv(fly,goodTracks(track),startNdx(i):endNdx(i))));
                    tmpCurve = mod(tmpCurve,360);
                elseif strcmpi(self.states.key{state},'curved walks')
                    tmpCurve = mean((allTracks.curv(fly,goodTracks(track),startNdx(i):endNdx(i)))).*self.fs;
                elseif strcmpi(self.states.key{state},'stops')
                    tmpCurve = abs(sum(allTracks.curv(fly,goodTracks(track),startNdx(i):endNdx(i))));
                    tmpCurve = mod(tmpCurve,360);
                elseif strcmpi(self.states.key{state},'boundary')
                    %tmpCurve = diff(abs(allTracks.curv(fly,goodTracks(track),[startNdx(i),endNdx(i)])));
                    tmpCurve = sum(allTracks.curv(fly,goodTracks(track),startNdx(i):endNdx(i)));
                end
                tmpSpd = mean(allTracks.spd(fly,goodTracks(track),startNdx(i):endNdx(i)));
                tmpDur = (endNdx(i)-startNdx(i)+1)./self.fs;
                
                if any(i == ndx)
                    stateDist{state}.curv = [stateDist{state}.curv, abs(tmpCurve)];
                    stateDist{state}.spd = [stateDist{state}.spd, tmpSpd];
                    stateDist{state}.dur = [stateDist{state}.dur, tmpDur];
                end
            end
        end
    end
end

end