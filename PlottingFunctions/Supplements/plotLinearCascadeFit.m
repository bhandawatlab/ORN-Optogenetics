function [fNum] = plotLinearCascadeFit(fNum,meta)
fs = 100;

% Intensity to LFP
%load([pwd '\Data\LFP_2_25.mat'],'LFP','LFP_fit','stim','b','optFitN','eta','rho','lamb')
load(string(meta.LFPFilterFile),'LFP','LFP_fit','stim','b','optFitN','eta','rho','lamb');
nStimPt = size(stim,2);
tt = (1:nStimPt)./fs;
tt_filt = (0:size(b,1)-2)./fs;

% define colors so that orange is left and blue is right
defColor = [0 0 0; 0 0.4470 0.7410];
%defColor = [0.8500 0.3250 0.0980; 0 0.4470 0.7410];

%% plot stim to LFP
figure(fNum);set(gcf,'Position',[2 42 838 924])
set(gcf,'defaultAxesColorOrder',defColor);
for i = 1:6
    subplot(6,1,i);
    yyaxis left; plot(tt,LFP(:,i),'-k','Linewidth',1);hold on;
    plot(tt,LFP_fit(:,i),'-r','Linewidth',1);
    ylabel('LFP (mv)');
    yyaxis right; plot(tt,stim(i,:),'Linewidth',1)
    xlabel('time (s)');ylabel('Stimulus (V)')
    
    if i == 1
        yyaxis left;
        legend({'empirical','model'})
    end
end

%% plot stim to LFP filter
figure(fNum+2);set(gcf,'Position',[2 42 838 924])
subplot(2,2,1);
plot(tt_filt,flip(b(2:end,optFitN)),'k','Linewidth',2);hold on;
plot([tt_filt(1) tt_filt(end)],[0 0],'k','Linewidth',1);
ylim([-0.01 0.04])
xlabel('lag (s)');title('LFP Filter')
subplot(2,2,2);
loglog(rho,eta);hold on;
plot(rho(optFitN),eta(optFitN),'or')
xlim([50 500]);ylim([10^-3 1])
text(100,0.01,['lambda = ' num2str(lamb(optFitN))])
text(100,0.01./2,['err = ' num2str(rho(optFitN))])
xlabel('residual norm ||Ax-b||_2');ylabel('solution norm ||x||_2')

%% plot LFP to firing rate
%load([pwd '\Data\LFP2Rate_4_12.mat'],'spkrate','spkrate_fit','stim','b','optFitN','eta','rho','lamb')
load(string(meta.RateFilterFile),'spkrate','spkrate_fit','stim','b','optFitN','eta','rho','lamb');
spkrate = reshape(spkrate,[],6)';

figure(fNum+1);set(gcf,'Position',[2 42 838 924])
set(gcf,'defaultAxesColorOrder',defColor);
for i = 1:6
    subplot(6,1,i);
    yyaxis left; plot(tt,spkrate(i,:),'-k','Linewidth',1);hold on;
    plot(tt,spkrate_fit(i,:),'-r','Linewidth',1);
    ylabel('firing rate (spikes/s)');
    yyaxis right; plot(tt,stim(i,:),'Linewidth',1)
    xlabel('time (s)');ylabel('Stimulus (V)')
    
    if i == 1
        yyaxis left;
        legend({'empirical','model'})
    end
end

%% plot LFP to firing rate filter
figure(fNum+2);set(gcf,'Position',[2 42 838 924])
subplot(2,2,3);
plot(tt_filt,flip(b(2:end,optFitN)),'k','Linewidth',2);hold on;
plot([tt_filt(1) tt_filt(end)],[0 0],'k','Linewidth',1);
ylim([-100 250]);xlim([0 0.5])
xlabel('lag (s)');title('rate Filter')
subplot(2,2,4);
loglog(rho,eta);hold on;
plot(rho(optFitN),eta(optFitN),'or')
xlim([700 2000]);ylim([1 1000])
text(1200,100,['lambda = ' num2str(lamb(optFitN))])
text(1200,100./2,['err = ' num2str(rho(optFitN))])
xlabel('residual norm ||Ax-b||_2');ylabel('solution norm ||x||_2')

fNum = fNum+3;

end