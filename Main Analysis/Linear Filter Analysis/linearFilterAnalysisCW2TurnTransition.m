function [b_leave,b_enter] = linearFilterAnalysisCW2TurnTransition(f_orco,transitionsOnly,meta,plotFigure)
% linearFilterAnalysisSpeedCurvature  wrapper function for calculating the
%   linear filters and turn triggered average for sharp turn transition
%   probability and probability of being in a turn at crossing (defined by
%   delta firing rate)
%
%   Inputs: f_orco = fly object
%           transitionsOnly = true/false (true = transition prob, false =
%               prob of being in sharp turn)
%           meta = structure with property border for light border (in cm)
%           plotFigure = true/false for plotting analysis figures
%
%   Output: b_leave = n x m linear filters for turn prob when leaving
%           b_enter = linear filters for turn prob when entering
%
close all
if plotFigure
    if ~exist(strcat(string(meta.plotFold),'/LinearFilterAnalysis/'), 'dir')
        mkdir(strcat(string(meta.plotFold),'/LinearFilterAnalysis/'))
    end
end
meta.delay = 10*30;
meta.buffer = 1*30;
meta.stNdx = find(strcmpi(f_orco.states.key  ,'sharp turns'));

%% turn triggered average
[turnNdxAbs,~,~] = turnTriggeredAverageAnalysis(f_orco,meta,plotFigure);

%% set up data for calculating average filter for all turns
%transitionsOnly = false;
%transitionsOnly = true;

if transitionsOnly
    yl = 'Turn transition Probability';
    figureFile = 'Turn transition Probability';
else
    yl = 'Probability in Turn State';
    figureFile = 'Probability in Turn State';
end

thresh = 15;

[df_leaving_all,df_entering_all,f_leaving_all,f_entering_all,turnProb_leaving,...
    turnProb_entering,turnProb_leaving_std,turnProb_entering_std] = ...
    getDf_aligned_turnRates(f_orco,thresh,meta.delay,turnNdxAbs,transitionsOnly,yl,f_orco.fs, plotFigure);

if plotFigure
    for i = 1:6
        subplot(3,2,i);
        xlim([-10 10])
    end
end

fs = 100;
newDelay = floor(meta.delay.*fs./f_orco.fs);
turnProb_leaving = resample(padarray(turnProb_leaving,[0 5],'replicate'),fs,f_orco.fs);
turnProb_leaving = turnProb_leaving((ceil(end./2)-newDelay):(ceil(end./2)+newDelay));
turnProb_entering = resample(padarray(turnProb_entering,[0 5],'replicate'),fs,f_orco.fs);
turnProb_entering = turnProb_entering((ceil(end./2)-newDelay):(ceil(end./2)+newDelay));
turnProb_leaving_std = resample(padarray(turnProb_leaving_std,[0 5],'replicate'),fs,f_orco.fs);
turnProb_leaving_std = turnProb_leaving_std((ceil(end./2)-newDelay):(ceil(end./2)+newDelay));
turnProb_entering_std = resample(padarray(turnProb_entering_std,[0 5],'replicate'),fs,f_orco.fs);
turnProb_entering_std = turnProb_entering_std((ceil(end./2)-newDelay):(ceil(end./2)+newDelay));

df_leaving_all = resample(padarray(df_leaving_all,[0 5],'replicate'),fs,f_orco.fs);
df_leaving_all = df_leaving_all((ceil(end./2)-newDelay):(ceil(end./2)+newDelay));
df_entering_all = resample(padarray(df_entering_all,[0 5],'replicate'),fs,f_orco.fs);
df_entering_all = df_entering_all((ceil(end./2)-newDelay):(ceil(end./2)+newDelay));
f_leaving_all = resample(padarray(f_leaving_all,[0 5],'replicate'),fs,f_orco.fs);
f_leaving_all = f_leaving_all((ceil(end./2)-newDelay):(ceil(end./2)+newDelay));
f_entering_all = resample(padarray(f_entering_all,[0 5],'replicate'),fs,f_orco.fs);
f_entering_all = f_entering_all((ceil(end./2)-newDelay):(ceil(end./2)+newDelay));


factors_leave_all = {df_leaving_all',df_entering_all',f_leaving_all',...
    f_entering_all',[df_leaving_all',f_leaving_all'],[df_entering_all',f_entering_all']};
turn_prob_all = {turnProb_leaving',turnProb_entering',turnProb_leaving',...
    turnProb_entering',turnProb_leaving',turnProb_entering'};
turn_prob_std_all = {turnProb_leaving_std',turnProb_entering_std',turnProb_leaving_std',...
    turnProb_entering_std',turnProb_leaving_std',turnProb_entering_std'};
title_id = {'Leaving - df fit','Entering - df fit','Leaving - f fit',...
    'Entering - f fit','Leaving - f+df fit','Entering - f+df fit'};

%% GLM analysis
if plotFigure
    figure;set(gcf,'Position',[2 42 838 924]);
    for i = 1:numel(factors_leave_all)
        subplot(3,2,i);
        [~,yhat(i,:)] = getLogisticRegression(factors_leave_all{i},turn_prob_all{i},turn_prob_std_all{i},newDelay,fs);
        xlabel('meta.delay');ylabel(yl);
        title( title_id{i})
    end
end

%% generate linear filters
% leaving
turnProb = turnProb_leaving(newDelay+1:end)';
turnProb_std = turnProb_leaving_std(newDelay+1:end)';
ntfilt = 2*fs;
nthist = 2*fs;

% firing rate filter
stim_all = {f_leaving_all};filter_type = {'firing rate filter'};
[U,s,V,XStimNew] = getSVD(stim_all,ntfilt,nthist,newDelay);
b_leave = generateLinearFilter_KinTurnProb(stim_all,turnProb,turnProb_std,ntfilt,filter_type,newDelay,yl,14,U,s,V,XStimNew,fs,plotFigure);
if plotFigure
    sgtitle('Leaving')
end

% entering
turnProb = turnProb_entering(newDelay+1:end)';
turnProb_std = turnProb_entering_std(newDelay+1:end)';
ntfilt = 2*fs;
nthist = 2*fs;

% firing rate filter
stim_all = {f_entering_all};filter_type = {'firing rate filter'};
[U,s,V,XStimNew] = getSVD(stim_all,ntfilt,nthist,newDelay);
b_enter = generateLinearFilter_KinTurnProb(stim_all,turnProb,turnProb_std,ntfilt,filter_type,newDelay,yl,14,U,s,V,XStimNew,fs,plotFigure);
if plotFigure
    sgtitle('Entering')
end

if plotFigure
    psFileName = strcat(string(meta.plotFold),'/LinearFilterAnalysis/', figureFile,'.ps');
    if exist(psFileName, 'file')==2
        delete(psFileName);
    end
    for f = 1:get(gcf,'Number')
        figure(f);
        print('-painters','-dpsc2',psFileName,'-loose','-append');
    end
    try
        ps2pdf('psfile', psFileName, 'pdffile', ...
            strcat(string(meta.plotFold),'/LinearFilterAnalysis/', figureFile,'.pdf'), 'gspapersize', 'letter',...
            'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
            'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
            'gslibpath','C:\Program Files\gs\gs9.50\lib');
    catch
        disp('No ghostscript available. Please install ghostscript or ')
        disp('change path to ghostscript in line 134 of linearFilterAnalysisCW2TurnTransition.m')
    end
end

end

function [turnNdxAbs,turnNdxAbs_Trans,turnNdxAbs_during] = turnTriggeredAverageAnalysis(f_orco,meta,plotFigure)
% compute turn triggered averages
border = meta.border;
df = f_orco.calcDeltaFR;
fe = f_orco.getFirstEntry('H',border);
fr_baseline = f_orco.spk(1);
stNdx = meta.stNdx;

turnNdxAbs = cell(1,f_orco.nFly);% ndx of peak in curvature during turn
turnNdxAbs_Trans = cell(1,f_orco.nFly);% ndx of transition into turn
turnNdxAbs_during = cell(1,f_orco.nFly);% entire turn
during_firingRate = [];
during_deltaFiringRate = [];
beforeDuring = false(f_orco.nFly,f_orco.nPt);
for fly = 1:f_orco.nFly
    [startNdx,endNdx,type] = startEndSeq(f_orco.states.ndx(fly,:)==stNdx);
    startNdx = startNdx(type==1);
    endNdx = endNdx(type==1);
    turnNdx = nan(1,numel(startNdx));
    maxCurv = nan(1,numel(startNdx));
    totCuv = nan(1,numel(startNdx));
    for track = 1:numel(startNdx)
        tmpCurv = f_orco.curv(fly,startNdx(track):endNdx(track));
        [maxCurv(track),turnNdx(track)] = max(tmpCurv);
        totCuv(track) = sum(tmpCurv).*180./pi;
    end
    turnNdxAbs{fly} = startNdx+turnNdx(track)-1;
    turnNdxAbs_Trans{fly} = startNdx;
    turnNdxAbs_during{fly} = find(f_orco.states.ndx(fly,:)==stNdx);


    during_firingRate = [during_firingRate f_orco.spk(fly,fe(fly):end)];
    during_deltaFiringRate = [during_deltaFiringRate df(fly,fe(fly):end)];
    beforeDuring(fly,fe(fly):end) = true;
end
meanFR_during = nanmean(during_firingRate);
meandFR_during = nanmean(during_deltaFiringRate);
meanFR_during_noBaseline = nanmean(during_firingRate(abs(during_firingRate-fr_baseline)>0.001));
meandFR_during_noBaseline = nanmean(during_deltaFiringRate(abs(during_firingRate-fr_baseline)>0.001));

turnNdx_cond = [turnNdxAbs;turnNdxAbs_during];
nCond = size(turnNdx_cond,1);
staByFly_f = cell(f_orco.nFly,nCond);
staByFly_df = cell(f_orco.nFly,nCond);
staByFly_beforeDuring = cell(f_orco.nFly,3);
staAllFly_f = nan(meta.delay,nCond);staAllFly_f_noBefore = nan(meta.delay,nCond);
staAllFly_df = nan(meta.delay,nCond);staAllFly_df_noBefore = nan(meta.delay,nCond);
staAllFly_f_prevEntry = nan(meta.delay,nCond); staAllFly_df_prevEntry = nan(meta.delay,nCond);
for cond = 1:nCond
    for fly = 1:f_orco.nFly
        currTurnNdxDuring = turnNdx_cond{cond,fly}(turnNdx_cond{cond,fly}>5400);
        curr_beforeDuring = beforeDuring(fly,:);
        if ~isempty(currTurnNdxDuring)
            curr_f = f_orco.spk(fly,:);
            curr_df = df(fly,:);
            staByFly_f{fly,cond} = curr_f(currTurnNdxDuring'+[-meta.delay:-1]);
            staByFly_df{fly,cond} = curr_df(currTurnNdxDuring'+[-meta.delay:-1]);
            staByFly_beforeDuring{fly,cond} = curr_beforeDuring(currTurnNdxDuring'+[-meta.delay:-1]);
        else
            staByFly_f{fly,cond} = nan(1,meta.delay);
            staByFly_df{fly,cond} = nan(1,meta.delay);
            staByFly_beforeDuring{fly,cond} = false(1,meta.delay);
        end
    end
    during = cell2mat(staByFly_beforeDuring(:,cond));
    sta_allF = cell2mat(staByFly_f(:,cond));
    staAllFly_f(:,cond) = nanmean(sta_allF);
    sta_allF(~during) = nan;
    staAllFly_f_noBefore(:,cond) = nanmean(sta_allF);

    sta_alldF = cell2mat(staByFly_df(:,cond));
    staAllFly_df(:,cond) = nanmean(sta_alldF);
    sta_alldF(~during) = nan;
    staAllFly_df_noBefore(:,cond) = nanmean(sta_alldF);

    for track = 1:size(sta_allF,1)
        ndx = find(abs(sta_allF(track,:)-fr_baseline)>0.001,1,'last');
        if isempty(ndx)
            ndx = 1;
        end
        sta_allF(track,ndx:end) = nan;
        sta_alldF(track,ndx:end) = nan;
    end
    staAllFly_f_prevEntry(:,cond) = nanmean(sta_allF);
    staAllFly_df_prevEntry(:,cond) = nanmean(sta_alldF);
end

if plotFigure
    figure;set(gcf,'Position',[2 42 838 924]);
    subplot(4,2,1);histogram(during_firingRate,'Normalization','probability');hold on;
    plot(meanFR_during.*[1 1],[0 0.3],'-k','LineWidth',2);
    plot(meanFR_during_noBaseline.*[1 1],[0 0.3],'--k','LineWidth',2);
    legend({'distribution','mean','mean no baseline'})
    xlabel('firing rate (spikes/s)');ylabel('Probability');
    title('firing rates after first entry')
    subplot(4,2,2);histogram(during_deltaFiringRate,'Normalization','probability');hold on;
    plot(meandFR_during.*[1 1],[0 0.3],'-k','LineWidth',2);xlim([-150 150])
    plot(meandFR_during_noBaseline.*[1 1],[0 0.3],'--k','LineWidth',2);xlim([-150 150])
    legend({'distribution','mean'})
    xlabel('delta firing rate (spikes/s^2)');ylabel('Probability');
    title('delta firing rates after first entry')
    subplot(4,2,3);plot([-meta.delay:-1]./30,staAllFly_f,'LineWidth',2);hold on;
    plot([-meta.delay 0]./30,meanFR_during.*[1 1],'-k','LineWidth',2);
    plot([-meta.delay 0]./30,meanFR_during_noBaseline.*[1 1],'--k','LineWidth',2);
    legend({'STA ST Curv Peak','STA ST duration','mean during','mean no baseline'},'Location','best')
    %legend({'STA ST Curv Peak','STA ST entry','STA ST duration','mean during','mean no baseline'},'Location','best')
    xlabel('meta.delay (s)');ylabel('firing rate (spikes/s)');
    title('STA includes before first entry (start of turns after first entry)')
    subplot(4,2,4);plot([-meta.delay:-1]./30,staAllFly_df,'LineWidth',2);hold on;
    plot([-meta.delay 0]./30,meandFR_during.*[1 1],'-k','LineWidth',2);
    plot([-meta.delay 0]./30,meandFR_during_noBaseline.*[1 1],'--k','LineWidth',2);
    xlabel('meta.delay (s)');ylabel('delta firing rate (spikes/s^2)');
    subplot(4,2,5);plot([-meta.delay:-1]./30,staAllFly_f_noBefore,'LineWidth',2);hold on;
    plot([-meta.delay 0]./30,meanFR_during.*[1 1],'-k','LineWidth',2);
    plot([-meta.delay 0]./30,meanFR_during_noBaseline.*[1 1],'--k','LineWidth',2);
    xlabel('meta.delay (s)');ylabel('firing rate (spikes/s)');
    title('STA excluding before first entry')
    subplot(4,2,6);plot([-meta.delay:-1]./30,staAllFly_df_noBefore,'LineWidth',2);hold on;
    plot([-meta.delay 0]./30,meandFR_during.*[1 1],'-k','LineWidth',2);
    plot([-meta.delay 0]./30,meandFR_during_noBaseline.*[1 1],'--k','LineWidth',2);
    xlabel('meta.delay (s)');ylabel('delta firing rate (spikes/s^2)');
    subplot(4,2,7);plot([-meta.delay:-1]./30,staAllFly_f_prevEntry,'LineWidth',2);hold on;
    plot([-meta.delay 0]./30,meanFR_during.*[1 1],'-k','LineWidth',2);
    plot([-meta.delay 0]./30,meanFR_during_noBaseline.*[1 1],'--k','LineWidth',2);
    xlabel('meta.delay (s)');ylabel('firing rate (spikes/s)');
    title('STA only considering until most recent entry')
    subplot(4,2,8);plot([-meta.delay:-1]./30,staAllFly_df_prevEntry,'LineWidth',2);hold on;
    plot([-meta.delay 0]./30,meandFR_during.*[1 1],'-k','LineWidth',2);
    plot([-meta.delay 0]./30,meandFR_during_noBaseline.*[1 1],'--k','LineWidth',2);
    xlabel('meta.delay (s)');ylabel('delta firing rate (spikes/s^2)');


    %print('-painters','-dpdf',strcat(string(meta.plotFold),'/LinearFilterAnalysis/STA_sharpTurnAnalysis.pdf'));
end
end

function [U,s,V,XStimNew] = getSVD(stim_all,ntfilt,nthist,delay)
opts.ncells = 1;
XStimNew_all = [];
for i = 1:numel(stim_all)
    stim = stim_all{i};
    % design matrix
    [XStimNew,~,ntfilt,nthist] = createDesignMat(stim,[],ntfilt,nthist,opts);
    if i == 1
        XStimNew_all = XStimNew;
    else
        XStimNew_all = [XStimNew_all,XStimNew(:,2:end)];
    end
end
XStimNew = XStimNew_all(delay+1:end,:);

% prep for the ridge regression using svd
tic;[U,s,V] = csvd(XStimNew,[]);toc
end