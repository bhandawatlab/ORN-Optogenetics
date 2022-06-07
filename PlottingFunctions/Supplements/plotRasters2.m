function [fNum] = plotRasters2(fNum,dataFile,functionFile,dataTraceFolder)
load(dataFile,'emp')
spikes = emp.rast;
fs = emp.fs1;
load(functionFile,'convVI','p')
I = convVI(emp.I);

nTrials_per_stim = size(spikes,2);
nTrials_tot = numel(spikes);
nStim = size(spikes,1);
tt = (1:size(I,2))./fs;

for i = 1:nTrials_tot
    spikes{i} = spikes{i}'./emp.fs2;
end

% plot supplementary figures
for i = 1:nStim
    s = dir([dataTraceFolder 'Type ' num2str(i) '/*.mat']);
    
    figure(fNum);set(gcf,'Position',[2 42 838 924])
    for j = 1:3
        load([dataTraceFolder 'Type ' num2str(i) '/' s(j).name]);
        subplot(4,1,j);hold on;
        plot((1:numel(DataExp.voltageBS))./emp.fs2,DataExp.voltageBS,'k');%ylim([-5 2])
        scatter(spikes{i,j},DataExp.voltageBS(round(spikes{i,j}.*emp.fs2)));hold off
        xlabel('time (s)');ylabel('Voltage (mV)')
        title(s(j).name(1:end-4), 'Interpreter', 'none')
        
        tmp = split(s(j).name(1:end-4),'_');
        trialN(i,j) = str2num(tmp{end});
    end
    subplot(4,1,4);
    plotSpikeRaster(spikes(i,:)','PlotType','vertline');hold on;
    plot(tt,-I(i,:),'k');ylim([-4 4])
    xlabel('time (s)');ylabel('Trials (Stimulus mW/cm^2)')
    fNum = fNum+1;
end

spikeType{1} = 1:6;
spikeType{2} = 7:12;
spikeType{3} = 13:18;
spikeTypeLabels  = {'ab2','unknown','ab1a'};

xx = [0:0.3:60];
for i = 1:nStim
    psth{i} = zeros(3,numel(xx));
    [trialN_sorted,ndx] = sort(trialN(i,:));
    for j = 1:3
        psth{i}(ndx(j),:) = histc(spikes{i,j},xx)./0.3;
    end
end

figure(fNum);set(gcf,'Position',[2 42 838 924])
cc = 'rck';
for i = 1:nStim
    %[~,ndx] = sort(trialN(i,:));
    subplot(6,1,i);hold on;
    for j = 1:3
        plot(xx,psth{i}(j,:),cc(j),'Linewidth',1);
    end
    xlim([0 60]);ylim([0 100])
	ylabel('# of spikes');xlabel('time (s)')
    legend(spikeTypeLabels)
    hold off;
end
fNum = fNum+1;

figure(fNum);set(gcf,'Position',[2 42 838 924])
figure(fNum+1);set(gcf,'Position',[2 42 838 924])
cc = 'rck';kk = 1;
correlation = zeros(nStim,3);
for i = 1:nStim
    [~,ndx] = sort(trialN(i,:));
    k = 1;
    for j = 1:2
        for jj = j+1:3
            tmp1 = psth{i}(j,:);
            tmp2 = psth{i}(jj,:);
            
            correlation(kk) = corr(tmp1',tmp2');
            RMSE = rms([tmp1-tmp2]);
            
            figure(fNum+1);
            subplot(3,1,k);hold on;
            scatter(tmp1,tmp2,'r');
            plot([0:100],[0:100],'--k','Linewidth',2)
            axis([0 100 0 100]);axis square
            set(gca, 'XTick', 0:25:100, 'YTick', 0:25:100);
            ylabel([spikeTypeLabels{jj} ' # of spikes']);
            xlabel([spikeTypeLabels{j} ' # of spikes'])
            k = k+1;
            axis square
            
            figure(fNum);
            subplot(6,3,kk);hold on;
            scatter(tmp1,tmp2,'r');
            plot([0:100],[0:100],'--k','Linewidth',2)
            axis([0 100 0 100]);axis square
            set(gca, 'XTick', 0:25:100, 'YTick', 0:25:100);
            ylabel([spikeTypeLabels{jj} ' # of spikes']);
            xlabel([spikeTypeLabels{j} ' # of spikes'])
            title(['Corr=' num2str(round(correlation(kk),2)) ', R_sq=' num2str(round(correlation(kk).^2,2)) ', RMSE=' num2str(round(RMSE,1))])
            kk = kk+1;
            axis square
            hold off;
        end
    end
end
fNum = fNum+2;

xScatter = (rand(6,3)-0.5)./2+repmat([1:3],6,1);
%X = categorical({'ab2, unknown','ab1a, ab2','ab1a, unknown'});
figure(fNum);
bar(1:3,mean(correlation));hold on
scatter(xScatter(:),correlation(:),'k')
er = errorbar(1:3,mean(correlation),std(correlation),std(correlation));
er.Color = [0 0 0];
er.LineStyle = 'none';
xticks([1:3])
xticklabels({'ab2, unknown','ab1a, ab2','ab1a, unknown'})
fNum = fNum+1;


end