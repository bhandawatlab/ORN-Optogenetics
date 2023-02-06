function [fNum] = plotTime2Return(f_orcoAll,border,fNum)

dur_baseline = cell(1,numel(f_orcoAll));
for i = 1:numel(f_orcoAll)
    fs = f_orcoAll{i}.fs;
    spk = f_orcoAll{i}.spk;
    dSpk = f_orcoAll{i}.calcDeltaFR;
    baseline = mean(spk(:,30*fs:60*fs),'all');
    fe = f_orcoAll{i}.getFirstEntry('H',border);
    
    % separate out into each time point into before first entry (FE), below
    % baseline firing rate, baseline firing rate after FE, and above
    % baseline firing rate
    condNdx = zeros(size(spk));condNdx(spk>baseline) = 4;
    condNdx(spk<baseline) = 2;condNdx(abs((spk-baseline))<0.001 & dSpk==0) = 3;
    dur_baseline{i} = [];
    for j = 1:f_orcoAll{i}.nFly
        condNdx(j,1:fe(j)-1) = 1;
        [startNdx,endNdx,type] = startEndSeq(condNdx(j,:)==3);
        startNdx = startNdx(type);
        endNdx = endNdx(type);
        
        dur_baseline{i} = [dur_baseline{i},endNdx-startNdx+1];
    end
end

opts.plotStaggered = true;
for i = 1:numel(f_orcoAll)
    if mod(i,12)==1
        figure(fNum);set(gcf,'Position',[2 42 838 924])
        sgtitle({'Time spent at baseline/exit';
            'thresholds: |spd-baseline|<0.001, |dSpk|=0'})
        k = 1;fNum = fNum+1;
    end
    opts.yl = [0 ceil(max(cell2mat(dur_baseline))./fs./20).*20];
    
    if k == 1
        plotLegend = true;
    else
        plotLegend = false;
    end
    subplot(4,3,k);
    violinPlotsStats(dur_baseline{i}'./fs,opts);
    %violinPlots(dur_baseline{i}./fs,[0 yl],true,plotLegend);
    title(f_orcoAll{i}.id)
    ylabel('time (seconds)');
    k = k+1;
end
end