function [] = plotTime2FirstStateTransition(f_orco,plotFig)
% calculate smoothed change in firing rate
dF = f_orco.calcDeltaFR;%gradient(f_orco.spk);
dF(:,end-1:end) = repmat(dF(:,end-2),1,2);
thresh =50;

crossingtype = 'exit';
[~,transitionProb_exit,transProbRand_exit,rast_exit,~] = getExitTrajectory(f_orco,crossingtype,dF,thresh);

crossingtype = 'enter';
[~,transitionProb_enter,transProbRand_enter,rast_enter,~] = getExitTrajectory(f_orco,crossingtype,dF,thresh);
transitionProb = {transitionProb_exit,transitionProb_enter};
transProbRand = {transProbRand_exit,transProbRand_enter};

crossingTypes = {'exit','enter'};
states = f_orco.states.key(1:2);
k = 1;
if plotFig
    figure;set(gcf,'Position',[2 42 838 924])
    for trans = 1:2
        for state = 1:2
            subplot(3,2,k);
            plot((0:1:90)./f_orco.fs, cumsum(transitionProb{trans}(state,:)),'k');hold on;
            plot((0:1:90)./f_orco.fs, cumsum(transProbRand{trans}(state,:)),'g');
            ylim([0 1])
            xlabel('time since peak (s)');ylabel('% of crossing tracks transitioning');
            legend({'During emp','Before Rand'})
            title([crossingTypes{trans} ' ' states{state} ' at first transition'])
            k = k+1;
        end
    end
end
subplot(3,2,5);
imagesc((1:size(rast_exit{1},2))./f_orco.fs,1:size(rast_exit{1},1),1-rast_exit{1});colormap(gray)
xlabel('time since negative peak');ylabel('peak #');
subplot(3,2,6);
imagesc((1:size(rast_enter{1},2))./f_orco.fs,1:size(rast_enter{1},1),1-rast_enter{1});colormap(gray)
xlabel('time since positive peak');ylabel('peak #');
end

function [exitTrajectory,transitionProb,transProbRand,rast,rastBaseline] = getExitTrajectory(f_orco,crossingtype,dF,thresh)

exitTrajectory = cell(f_orco.nFly,1);
for fly = 1:f_orco.nFly
    if strcmpi(crossingtype,'enter')
        % positive dF
        [pks,locs] = findpeaks([dF(fly,:),0],'MinPeakHeight',thresh,'MinPeakDistance',10.*f_orco.fs./30);
        % negative dF
        [pks_p,locs_p] = findpeaks(-[dF(fly,:),0],'MinPeakHeight',thresh,'MinPeakDistance',10.*f_orco.fs./30);
    elseif strcmpi(crossingtype,'exit')
        % negative dF
        [pks,locs] = findpeaks(-[dF(fly,:),0],'MinPeakHeight',thresh,'MinPeakDistance',10.*f_orco.fs./30);
        % positive dF
        [pks_p,locs_p] = findpeaks([dF(fly,:),0],'MinPeakHeight',thresh,'MinPeakDistance',10.*f_orco.fs./30);
        
    end
    locs(locs>10700.*f_orco.nPt./10800) = [];
    locs_p(locs_p>10700.*f_orco.nPt./10800) = [];
    
    
    %     % plotting to check
    %     %--------------------
    %     figure(1);subplot(2,1,1);
    %     plot(dF(fly,:));hold on;
    %     scatter(locs,-pks,'r');scatter(locs_p,pks_p,'g');
    %     xlim([5400 10800]);hold off
    %     subplot(2,1,2);plot(f_orco.spk(fly,:));hold on;
    %     scatter(locs,f_orco.spk(fly,locs),'r');scatter(locs_p,f_orco.spk(fly,locs_p),'g');
    %     xlim([5400 10800]);hold off
    %     %--------------------
    
    % get exit trajectories
    [~,closestPtAft,idxB] = findBeforeAfter(locs,[locs,locs_p],'after');
    exitTrajectory{fly,1} = [locs;closestPtAft]';
end

% at crossing
allState = [];
for fly = 1:f_orco.nFly
    stateTmp = f_orco.states.ndx(fly,:);
    allState = [allState;stateTmp(exitTrajectory{fly,1}(:,1)+repmat([0:90],numel(exitTrajectory{fly,1}(:,1)),1))];
end

for i = 1:2
    allStateCross = allState(allState(:,1)==i,:);%61
    [~,timeSincePeak]=max(allStateCross~=i,[],2);
    timeSincePeak(timeSincePeak==1) = 91;
    transitionProb(i,:) = histcounts(timeSincePeak,[-0.5:1:90.5])./numel(timeSincePeak);

    rast{i} = zeros(numel(timeSincePeak),max([timeSincePeak; 3*f_orco.fs]));
    for j = 1:numel(timeSincePeak)
        rast{i}(j,timeSincePeak(j)) = 1;
    end
end

% random from before
allTrackRand = [];allStateRand = [];
for fly = 1:f_orco.nFly
    stateTmp = f_orco.states.ndx(fly,:);
    allStateRand = [allStateRand;stateTmp(ceil(5000.*rand(200,1))+61+repmat([0:90],200,1))];
end

for i = 1:2
    allStateCrossRand = allStateRand(allStateRand(:,1)==i,:);%61
    [~,timeSincePeak]=max(allStateCrossRand~=i,[],2);
    %timeSincePeak(timeSincePeak==1) = 91;
    timeSincePeak(timeSincePeak==1,:) = [];
    transProbRand(i,:) = histcounts(timeSincePeak,[-0.5:1:90.5])./numel(timeSincePeak);
    
    rastBaseline{i} = zeros(numel(timeSincePeak),max([timeSincePeak; 3*f_orco.fs]));
    for j = 1:numel(timeSincePeak)
        rastBaseline{i}(j,timeSincePeak(j)) = 1;
    end

end

end

