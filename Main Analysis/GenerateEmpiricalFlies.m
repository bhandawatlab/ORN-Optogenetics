function [] = GenerateEmpiricalFlies(genAll,meta)

border = meta.border;
fs = meta.fs;
rBound = meta.rBound;

for i = 1 :numel(genAll)
    tic;
    gen = genAll{i};
    % load data
    load([pwd '\Data\DataGen\' gen '_March2022.mat'],'Data','curvPks','curvWalks','stopCond','boundCond')
    % load spike prediction
    load([pwd '\Data\DataSpike\' gen '_SpkRate.mat'],'sps_pred2')
    
    % get the baseline 
    baseline = mean(sps_pred2(:,30*fs:60*fs),'all');
    % remove any artifacts from the linear filter
    sps_pred2(:,1:30*fs) = round(baseline,4);
    sps_pred2 = round(sps_pred2,4);
    sps_pred2 = sps_pred2(:,1:size(Data.x,2));
    
    % create flies object
    f_orco = Flies(gen,Data.x,Data.y,Data.xHead,Data.yHead,sps_pred2,fs,rBound,...
        Data.lightOn,curvPks,curvWalks,stopCond,boundCond,[]);
    
    % calculate first entry
    fe = f_orco.getFirstEntry('H',border);

    % remove flies that do not enter
    f_orco = f_orco.rmvData(isnan(fe));
    fe = f_orco.getFirstEntry('H',border);
    
    %%%%
    % for Or42a only
    if strcmpi(gen,'Or42a Retinal')
        fe(6) = nan;
        f_orco = f_orco.rmvData(isnan(fe));
        fe = f_orco.getFirstEntry('H',border);
    end
    %%%%
    
    % get the firing rate and change in firing rate
    spk = f_orco.spk;
    dSpk = f_orco.calcDeltaFR;
    
    % separate out into each time point into before first entry (FE), below
    % baseline firing rate, baseline firing rate after FE, and above
    % baseline firing rate
    condNdx = zeros(size(spk));condNdx(spk>baseline) = 4;
    condNdx(spk<baseline) = 2;condNdx(abs((spk-baseline))<0.001 & dSpk==0) = 3;
    for j = 1:f_orco.nFly
        condNdx(j,1:fe(j)-1) = 1;
    end
    key = {'before','below','baseline','above'};
    
    % get the kinematics and decision space of the flies
    GetKinematicModelParams(f_orco,condNdx,key,fe,meta,true);
    
    % intensity space of the arena
    load('Data/Intensity_space2.mat','xN','p','Intensity_spaceN','convIV','convVI')
    f_orco.model.IntensitySpace.x = xN;
    f_orco.model.IntensitySpace.I = Intensity_spaceN;
    f_orco.model.IntensitySpace.p = p;
    f_orco.model.IntensitySpace.convIV = convIV;% conversion from intensity to voltage
    f_orco.model.IntensitySpace.convIV = convVI;% conversion from voltage to intensity
    
    % when df is 0 and there is a nonbaseline firing rate (inhibition period)
    %crossingtype = 'enter';
    crossingtype = 'exit';%plotFig = false;
    [f_orco,~,~,~,~,~,~,~] = getDistInhibition(f_orco,crossingtype,meta.plotFig);
    
%     % get border choice (prob of turn) due to large changes in firing rate
%     state = 2;% curved walk
%     [f_orco] = getBorderTurnRate(f_orco,state,false);
%     state = 1;% sharp turn
%     [f_orco] = getBorderTurnRate(f_orco,state,false);
    
    if meta.plotFig
        fName = ['Figures\' gen '_' meta.d meta.ext '_DataGen'];
        fNum = get(gcf,'Number');
        for f = 1:fNum-1
            figure(f);
            print('-painters','-dpsc2',[fName '.ps'],'-loose','-append');
        end
        ps2pdf('psfile', [fName '.ps'], 'pdffile', [fName '.pdf'], 'gspapersize', 'letter',...
            'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
            'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
            'gslibpath','C:\Program Files\gs\gs9.50\lib');
    end
    
    save(['DataModel/' gen '_' meta.d meta.ext '.mat'],'f_orco');%_allTime
    close all
    fprintf('Finished generating %s object in %d seconds\n',gen,round(toc))
end

end














