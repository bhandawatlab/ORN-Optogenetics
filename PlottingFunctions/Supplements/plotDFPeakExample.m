function [fNum] = plotDFPeakExample(f_orco,fNum)

% calculate smoothed change in firing rate
dF = f_orco.calcDeltaFR;%gradient(f_orco.spk);
dF(:,end-1:end) = repmat(dF(:,end-2),1,2);

thresh =50;
% peaks in df
NegDF = [];PosDF = [];dt = 60;%2s
for fly = 1:f_orco.nFly
    % negative dF
    [pks_n,locs_n] = findpeaks(-[dF(fly,:),0],'MinPeakHeight',thresh,'MinPeakDistance',15.*f_orco.fs./30);
    % positive dF
    [pks_p,locs_p] = findpeaks([dF(fly,:),0],'MinPeakHeight',thresh,'MinPeakDistance',15.*f_orco.fs./30);
    
    % plot a sample track
    if fly == 1
        %--------------------
        figure(fNum);set(gcf,'Position',[2 42 838 924])
        subplot(2,1,1);
        plot([1:f_orco.nPt]./f_orco.fs, dF(fly,:),'k','Linewidth',1);hold on;
        scatter(locs_n./f_orco.fs,-pks_n,'r','Linewidth',1);
        scatter(locs_p./f_orco.fs,pks_p,'g','Linewidth',1);
        xlim([180 360]);ylim([-200 200]);hold off
        xlabel('time (s)');ylabel('dF (spks/s^2)')
        subplot(2,1,2);
        plot([1:f_orco.nPt]./f_orco.fs, f_orco.spk(fly,:),'k','Linewidth',1);hold on;
        scatter(locs_n./f_orco.fs,f_orco.spk(fly,locs_n),'r','Linewidth',1);
        scatter(locs_p./f_orco.fs,f_orco.spk(fly,locs_p),'g','Linewidth',1);
        xlim([180 360]);ylim([0 45]);hold off
        xlabel('time (s)');ylabel('F (spks/s)')
        %--------------------
        fNum = fNum+1;
    end
    
    locs_n(locs_n>(f_orco.nPt-dt)) = [];
    locs_p(locs_p>(f_orco.nPt-dt)) = [];
    dFTmp = dF(fly,:);
    NegDF = [NegDF;dFTmp(locs_n'+[-dt:dt])];
    PosDF = [PosDF;dFTmp(locs_p'+[-dt:dt])];
end

% get the half width of the average trace
midPt = dt+1;
avgNegDF = mean(NegDF);
idxaft = find(avgNegDF(midPt:end)>avgNegDF(midPt)./2,1,'first')+midPt;
idxbef = find(avgNegDF(1:midPt)>avgNegDF(midPt)./2,1,'last');
halfWidthNeg = (idxaft-idxbef+1)./30;

avgPosDF = mean(PosDF);
idxaft = find(avgPosDF(midPt:end)<avgPosDF(midPt)./2,1,'first')+midPt;
idxbef = find(avgPosDF(1:midPt)<avgPosDF(midPt)./2,1,'last');
halfWidthPos = (idxaft-idxbef+1)./30;

%% plotting functions
% plot the average peak and the standard deviation around the average
figure(fNum);set(gcf,'Position',[2 42 838 924])
subplot(2,1,1);
shadedErrorBar((-dt:dt)./f_orco.fs,mean(NegDF),std(NegDF),'lineprops','k')
ylim([-200 100])
xlabel('time since peak');ylabel('dF (spks/s^2)');
legend({'stand dev','','','mean'})
title(['average half-width= ' num2str(halfWidthNeg) ' s'])
subplot(2,1,2);
shadedErrorBar((-dt:dt)./f_orco.fs,mean(PosDF),std(PosDF),'lineprops','k')
ylim([-100 200])
xlabel('time since peak');ylabel('dF (spks/s^2)');
title(['average half-width= ' num2str(halfWidthPos) ' s'])

fNum = fNum+1;
end



