function [f_orco_all] = GenerateEmpiricalFlies(genAll,meta)
% GenerateEmpiricalFlies  Generates a fly object, then calculate the KNN
%   space, baseline and boundary conditions, and crossing turn rate/
%   kinematic filters
%
%   Inputs: genAll = cell array of genotype labels to create fly objects for
%           meta = a structure with meta information
%
%   Output: f_orco_all = cell array of fly objects, where each cell is a
%           different genotype
%   

border = meta.border;
fs = meta.fs;
rBound = meta.rBound;
plotFigure = false;
genAll(cellfun(@(x) isempty(x),genAll)) = [];

f_orco_all = cell(1,numel(genAll));
for i = 1:numel(genAll)
    tic;
    gen = genAll{i};
    % load data
    load(strcat(string(meta.folderData),'\',gen,'_',meta.d, '.mat'),'Data','curvPks','curvWalks','stopCond','boundCond')
    % load spike prediction
    load(strcat(string(meta.foldSpk),'\',gen,'_SpkRate.mat'),'sps_pred2','insideSS','baseline')
    
    % create flies object
    f_orco = Flies(gen,Data.x,Data.y,Data.xHead,Data.yHead,sps_pred2,...
        insideSS,baseline,fs,rBound,Data.lightOn,curvPks,curvWalks,...
        stopCond,boundCond,[]);
    
    spd = f_orco.calcSpd;
    f_orco.states.ndx(spd<meta.stopThresh) = 3;
    
    % calculate first entry
    fe = f_orco.getFirstEntry('H',border);

    % remove flies that do not enter
    f_orco = f_orco.rmvData(isnan(fe));
    fe = f_orco.getFirstEntry('H',border);
    
    % get the firing rate and change in firing rate
    spk = f_orco.spk;
    dSpk = f_orco.calcDeltaFR;
    
    % separate out into each time point into before first entry (FE), below
    % baseline firing rate, baseline firing rate after FE, and above
    % baseline firing rate
    condNdx = zeros(size(spk));condNdx(spk>baseline) = 4;
    condNdx(spk<baseline) = 2;condNdx(abs((spk-baseline))<0.01 & dSpk<0.1) = 3;
    for j = 1:f_orco.nFly
        condNdx(j,1:fe(j)-1) = 1;
    end
    key = {'before','below','baseline','above'};
    
    % get the kinematics and decision space of the flies
    GetKinematicModelParams(f_orco,condNdx,key,fe,meta,plotFigure);
    
    % intensity space of the arena
    load(strcat(meta.DestPath,'/Data/Intensity_space2.mat'),'xN','p','Intensity_spaceN','convIV','convVI')
    f_orco.model.IntensitySpace.x = xN;
    f_orco.model.IntensitySpace.I = Intensity_spaceN;
    f_orco.model.IntensitySpace.p = p;
    f_orco.model.IntensitySpace.convIV = convIV;% conversion from intensity to voltage
    f_orco.model.IntensitySpace.convIV = convVI;% conversion from voltage to intensity
    
    % when df is 0 and there is a nonbaseline firing rate (inhibition period)
    %crossingtype = 'enter';
    crossingtype = 'exit';
    [f_orco,~,~,~,~,~,~] = getDistInhibition(f_orco,crossingtype,meta.timeInterval,plotFigure);
    
    % speed, curvature, and angle when flies are leaving the arena boundary
    f_orco = getLeavingBoundaryParam(f_orco,fe,plotFigure);
    
    % probability of initiating a sharp turn linear filter when entering 
    % and leaving the light border
    transitionOnly = true;
    [b_leave,b_enter] = linearFilterAnalysisCW2TurnTransition(f_orco,transitionOnly,meta,plotFigure);
    f_orco.model.border.leave.ST_trans = b_leave(:,14);
    f_orco.model.border.enter.ST_trans = b_enter(:,14);
    % kinematics linear filters when entering and leaving the light border
    [b_spdLeave,b_curvLeave,b_spdEnter,b_curvEnter] = linearFilterAnalysisSpeedCurvature(f_orco,meta,plotFigure);
    f_orco.model.border.leave.speed = b_spdLeave(:,14);
    f_orco.model.border.leave.curve = b_curvLeave(:,14);
    f_orco.model.border.enter.speed = b_spdEnter(:,14);
    f_orco.model.border.enter.curve = b_curvEnter(:,14);
    
    % older analysis, not used
%     % get border choice (prob of turn) due to large changes in firing rate
%     state = 2;% curved walk
%     [f_orco] = getBorderTurnRate(f_orco,state,false);
%     state = 1;% sharp turn
%     [f_orco] = getBorderTurnRate(f_orco,state,false);
    
%     % printing to pdf file
%     if meta.plotFig
%         fName = strcat(string(meta.plotFold),gen,'_',meta.d,meta.ext,'_DataGen');
%         fNum = get(gcf,'Number');
%         for f = 1:fNum-1
%             figure(f);
%             print('-painters','-dpsc2',[fName '.ps'],'-loose','-append');
%         end
%         ps2pdf('psfile', [fName '.ps'], 'pdffile', [fName '.pdf'], 'gspapersize', 'letter',...
%             'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
%             'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
%             'gslibpath','C:\Program Files\gs\gs9.50\lib');
%     end
    
    % saving the fly object
    if meta.saveData
        save(strcat(string(meta.foldDataModel),'\',gen,'_',meta.d,meta.ext,'.mat'),'f_orco');%_allTime
    end
    close all
    fprintf('Finished generating %s object in %d seconds\n',gen,round(toc))
    
    f_orco_all{i} = f_orco;
end

end














