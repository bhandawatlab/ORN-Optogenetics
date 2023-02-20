function PlotAllFiguresEmp(genAll,genRetinal,meta)

plotFold = [char(meta.plotFold) '\'];
subFoldName = meta.foldName;

% load in all data (both retinal and controls)
%--------------------------------------------------------------------------
% (for both control and retinal)
f_orcoAll = cell(1,numel(genAll));
for i = 1:numel(genAll)
    gen = genAll{i};
    try
        %load(['DataModel/' gen '_' meta.d meta.ext '.mat'],'f_orco');
        load(strcat(string(meta.foldDataModel),'\',gen,'_',meta.d,meta.ext,'.mat'),'f_orco');
        f_orcoAll{i} = f_orco;
    catch
        f_orcoAll{i} = [];
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting analysis figures for retinal and controls %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% plot time 2 return
close all;fNum = 1;
fNum = plotTime2Return(f_orcoAll,meta.border,fNum);
printFigures(fNum,[plotFold 'GeneralBehavior/'],'Time2Return_allGenotypes')
%%
% plot probability of being inside
close all;fNum = 1;
for i = 1:numel(f_orcoAll)
    if mod(i,10)==1
        figure(fNum);set(gcf,'Position',[2 42 838 924])
        k = 1;fNum = fNum+1;
    end
    subplot(5,2,k);
    tmpFly.rHead = f_orcoAll{i}.rH;
    plottingProbInside(tmpFly,meta.border,f_orcoAll{i}.nFly,...
        f_orcoAll{i}.fs,f_orcoAll{i}.lightOn,[]);
    title(f_orcoAll{i}.id);
    xlabel('time (minutes)');
    ylabel('Probability inside')
    k = k+1;
end
printFigures(fNum,[plotFold 'GeneralBehavior/'],'ProbabilityInside_allGenotypes')

%%
% plot probability of being inside (since first entry)
close all;fNum = 1;
for i = 1:numel(f_orcoAll)
    if mod(i,10)==1
        figure(fNum);set(gcf,'Position',[2 42 838 924])
        k = 1;fNum = fNum+1;
    end
    subplot(5,2,k);
    tmpFly.rHead = f_orcoAll{i}.rH;
    fe = f_orcoAll{i}.getFirstEntry('H',meta.border);
    sliding_mean{i} = plottingProbInside(tmpFly,meta.border,f_orcoAll{i}.nFly,...
        f_orcoAll{i}.fs,fe,'during');
    title(f_orcoAll{i}.id);
    xlabel('time (minutes)');
    ylabel('Probability inside')
    k = k+1;
end
printFigures(fNum,[plotFold 'GeneralBehavior/'],'ProbabilityInside_allGenotypes_AfterFE')
%% 
% plot tracks
close all;
for i = 1:numel(f_orcoAll)
    plotXYTracks(f_orcoAll{i},meta.border,f_orcoAll{i}.nFly,true,'')
    printFigures(get(gcf,'Number')+1,[plotFold 'GeneralBehavior/XY tracks/'],[f_orcoAll{i}.id ' Tracks'])
    close all;
end
%% 
% plot radial occupancy (for both control and retinal)
close all;fNum = 1;
stopSpd = 0;%meta.stopThresh;
[fNum] = plotRadialOccupancy(f_orcoAll,'H',meta.border,stopSpd,0,fNum);% plots radial occupancy with stops
printFigures(fNum,[plotFold 'GeneralBehavior/'],'RadialOccupancy_allGenotypes')

% plot radial occupancy (for both control and retinal)
close all;fNum = 1;
stopSpd = meta.stopThresh;
[fNum] = plotRadialOccupancy(f_orcoAll,'H',meta.border,stopSpd,0,fNum);% plots radial occupancy with stops
printFigures(fNum,[plotFold 'GeneralBehavior/'],'RadialOccupancy_noStops_allGenotypes')

%%
% plot spatial-temporal position comparison between controls and retinal
close all;fNum = 1;
for i = 1:numel(f_orcoAll)
    if mod(i,6)==1
        figure(fNum);set(gcf,'Position',[2 42 1000 924])
        fNum = fNum+1;k = 1;
    end
    subplot(3,2,k);
    plotRadialPosition(f_orcoAll{i},meta.border,'after',true,'dt',...
        2.*f_orcoAll{i}.fs,'clims',[0 0.01]);
    title(f_orcoAll{i}.id)
    if mod(k,2)==0
        ax2 = gca;
        currCLim(2,:) = ax2.CLim;
        maxCLim = [0 max(currCLim(:,2))];
        ax1.CLim = maxCLim;
        ax2.CLim = maxCLim;
    else
        ax1 = gca;
        currCLim(1,:) = ax1.CLim;
    end
    k = k+1;
    
    if i == numel(f_orcoAll) || k==7
        printFigures([fNum-1 nan],[plotFold 'GeneralBehavior/'],...
            ['SpatialTemporal_Density_RetinalControl' num2str(fNum-1)])
    end
end

% plot spatial-temporal position comparison for retinal
[retinalNdx,~]=ismember(genAll,genRetinal);
f_orcoRetinal = f_orcoAll(retinalNdx);
close all;fNum = 1;
for i = 1:numel(f_orcoRetinal)
    if mod(i,4)==1
        figure(fNum);set(gcf,'Position',[2 42 1000 924])
        fNum = fNum+1;k = 1;
    end
    subplot(2,2,k);
    plotRadialPosition(f_orcoRetinal{i},meta.border,'after',true,'dt',...
        2.*f_orcoRetinal{i}.fs,'clims',[0 0.01]);
    title(f_orcoRetinal{i}.id)
    ax = gca;
    ax.CLim = [0 ax.CLim(end)];
    k = k+1;
    if i == numel(f_orcoRetinal) || k==5
        printFigures([fNum-1 nan],[plotFold 'GeneralBehavior/'],...
            ['SpatialTemporal_Density_Retinal' num2str(fNum-1)])
    end
end

%printFigures(fNum,[plotFold 'GeneralBehavior/'],'SpatialTemporal_Density_AllGenotypes')

%%
% plot KS test of retinal vs control (only for time averaged)
if meta.adaptation == false
    close all;fNum = 1;
    [fNum,~,~] = plotStatTestRetContKNN(reshape(f_orcoAll,2,[]),meta.States2Plot_KNN,0,0.05,fNum);
    printFigures(fNum,[plotFold meta.foldName],'KNN_Control_Retinal_KS-test_p_05_color')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting analysis figures for control only %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[retinalNdx,~]=ismember(genAll,genRetinal);
controlNdx = ~retinalNdx;

%%
% plot KNN kinematic distributions image and relative colormap for control flies
close all;fNum = 1;
[~,currState,currKin] = plotKNNKinematics2(f_orcoAll(controlNdx),meta.States2Plot_KNN,...
    meta.tSlice,'imagesc','relative',fNum);
plotKNNDistributions(currState,currKin,[plotFold subFoldName 'KNN_Relative_Heatmap/'],...
    '_relative_heatmap_allGenotypes_Control')

%%
% plot turn optimality imagesc and relative colormap for control flies
close all;fNum = 1;
[fNum,currState] = plotTurnOpt2(f_orcoAll(controlNdx),meta.States2Plot_Opt,meta.tSlice,...
    'imagesc','relative',true,false,fNum);
currKin = cellstr(repmat('TurnOptimality',numel(currState),1))';
plotKNNDistributions(currState,currKin,[plotFold subFoldName 'KNN_Relative_Heatmap/'],'_relative_heatmap_Control')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting analysis figures for retinal only %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (for retinal only)
f_orcoAll = cell(1,numel(genRetinal));
for i = 1:numel(genRetinal)
    gen = genRetinal{i};
    try
        %load(['DataModel/' gen '_' meta.d meta.ext '.mat'],'f_orco');
        load(strcat(string(meta.foldDataModel),'\',gen,'_',meta.d,meta.ext,'.mat'),'f_orco');
        % perform permutation test across time slices for retinal flies (only for adaptation)
        if meta.adaptation == true
            [f_orco] = calcKNN_Habituation_PermutationTest(f_orco,1:4);
        end
        f_orcoAll{i} = f_orco;
    catch
        f_orcoAll{i} = [];
    end
end

% close all;fNum = 1;
% [fNum] = compareKNN(f_orcoAll,fNum);
% printFigures(fNum,'Figures/','KNN_3_OR_smaller_v_larger_dataset')

% %%
% % create subfolders for plots
% KNN_folders = {'KNN_Absolute_Heatmap','KNN_Relative_Heatmap','Inhibition'};
% subFoldName = meta.foldName;
% for i = 1:numel(KNN_folders)
%     if ~exist(['Figures/' foldName KNN_folders{i}], 'dir')
%         mkdir(['Figures/' foldName KNN_folders{i}])
%     end
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting one time analysis figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if meta.plotSupplements == true
    %%
    % plot sample trajectory in f/df space and distribution of data in f/df space
    close all;fNum = 1;
    [fNum] = plotF_df_SpacePlots(f_orcoAll,fNum,meta.border);
    printFigures(fNum,[plotFold 'f_dfInfo/'],'F_dF_Dist_AllGenotypes2')
    
    %
    % plot Rasters and ORN fits
    close all;fNum = 1;
    [fNum] = plotRasters2(fNum,meta.spikingDataFile,meta.Intensity2VoltageFile,meta.trainingDataFolder);
    printFigures(fNum,[plotFold 'LinearFilterAnalysis/'],'Spike rasters Filter')
    
    %%
    % plot Linear-linear cascade fit
    close all;fNum = 1;
    [fNum] = plotLinearCascadeFit(fNum,meta);
    printFigures(fNum,[plotFold 'LinearFilterAnalysis/'],'Linear Linear Filter')
    
    %%
    % plot dF peak schematic
    close all;fNum = 1;
    for i = 1:numel(f_orcoAll)
        gen = f_orcoAll{i}.id;
        fNumOld = fNum;
        [fNum] = plotDFPeakExample(f_orcoAll{i},fNum);
        for f = fNumOld:(fNum-1)
            figure(f)
            suptitle(gen)
        end
    end
    printFigures(fNum,[plotFold 'f_dfInfo/'],'Peak_delta_firing_rate_allGenotypes')
    
    %%
    % plot KNN formulation
    close all;fNum = 1;
    for i = 1:numel(f_orcoAll)
        [fNum] = plotKNNFormulation(f_orcoAll{i},fNum);
        fNum = fNum-1;
    end
    printFigures(fNum+1,[plotFold 'f_dfInfo/'],'KNN_formulation_allGenotypes')
    
    %%
%     % plot time at steady state
%     close all;fNum = 1;
%     fNum = plotTimeSteadyState(f_orcoAll,meta.border,fNum);
%     printFigures(fNum,'Figures/','TimeAtSteadyStateVsCrossing_allGenotypes')
    
    %%
    % plot turn triggered average
    close all;fNum = 1;
    plotSchematicTTA(f_orcoAll{1});
    printFigures(fNum,[plotFold 'LinearFilterAnalysis/'],'TurnTriggeredAverageExampleTracks')
    
    %%
    % plot time to transition
    close all;fNum = 1;
    plotTime2FirstStateTransition(f_orcoAll{1},true);
    printFigures(fNum,[plotFold 'f_dfInfo/'],'Time2FirstStateTransitionAfterPeakInDf')
    
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting KNN based analysis figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% plot permutation test across time slices for retinal flies (only for adaptation)
if meta.adaptation == true
    close all;fNum = 1;
    %%% plot KS test across time slices for retinal flies (only for
    %%% adaptation) - do not use
    %%%[~,currState,currKin] = plotStatTestAdaptationKNN(f_orcoAll,meta.States2Plot_KNN,meta.tSlice,fNum);
    [~,currState,currKin] = plotKNN_Habituation_PermutationTest...
        (f_orcoAll,1:4,meta.tSlice,fNum);
    plotKNNDistributions(currState,currKin,[plotFold subFoldName 'KNN_Emp_KS_test/'],...
        '__allGenotypes')
end

%%
% plot comparison of different distribution fits
if meta.adaptation == false
    close all;fNum = 1;
    fNum = plotDistributionComparisons(f_orcoAll,meta.States2Plot_KNN,fNum);
    %[fNum,~,~] = plotDistributionComparisons2(f_orcoAll,meta.States2Plot_KNN,fNum);
    printFigures(fNum,[plotFold subFoldName],'KNN_Distribution_Fit_2')
end

%%
% plot transition probability image and relative colormap
close all;fNum = 1;
[fNum] = plotTransitionProbability2(f_orcoAll,'imagesc','relative',fNum);
printFigures(fNum,[plotFold subFoldName],'TransitionProbability_relative_heatmap_allGenotypes')

%%
% plot transition probability image and absolute colormap
close all;fNum = 1;
[fNum] = plotTransitionProbability2(f_orcoAll,'imagesc','absolute',fNum);
printFigures(fNum,[plotFold subFoldName],'TransitionProbability_absolute_heatmap_allGenotypes')

%%
% plot baseline kinematic comparisons
close all;fNum = 1;
[fNum] = plotBaselineKinematics(f_orcoAll,meta.border,fNum);
printFigures(fNum,[plotFold subFoldName],'BaselineKinematics')

%%
% plot spatial probability map
close all;fNum = 1;
[fNum] = plotROISpatialProb(f_orcoAll,[0.1 0.1],'H',meta.border,fNum);
printFigures(fNum,[plotFold subFoldName],'SpatialProbMap')

%
% plot inhibition kinematics
close all;fNum = 1;
[~,currState,currKin] = plotInhibitionKinematics2(f_orcoAll,...
    meta.States2Plot_Inh,meta.tSlice,fNum);
plotKNNDistributions(currState,currKin,[plotFold subFoldName 'Inhibition/'],...
    '_inhibition_allGenotypes')

%%
% plot KNN kinematic distributions image and absolute colormap
close all;fNum = 1;
[~,currState,currKin] = plotKNNKinematics2(f_orcoAll,meta.States2Plot_KNN,...
    meta.tSlice,'imagesc','absolute',fNum);
plotKNNDistributions(currState,currKin,[plotFold subFoldName 'KNN_Absolute_Heatmap/'],...
    '_absolute_heatmap_allGenotypes')

%%
% plot KNN kinematic distributions image and relative colormap
close all;fNum = 1;
[~,currState,currKin] = plotKNNKinematics2(f_orcoAll,meta.States2Plot_KNN,...
    meta.tSlice,'imagesc','relative',fNum);
plotKNNDistributions(currState,currKin,[plotFold subFoldName 'KNN_Relative_Heatmap/'],...
    '_relative_heatmap_allGenotypes')

%%
% plot KNN kinematic distributions standev image and absolute colormap
close all;fNum = 1;
[~,currState,currKin] = plotKNNKinematicsSTD(f_orcoAll,meta.States2Plot_KNN,...
    meta.tSlice,'imagesc','absolute',fNum);
plotKNNDistributions(currState,currKin,[plotFold subFoldName 'KNN_Absolute_Heatmap/'],...
    '_absolute_heatmap_allGenotypes_std')

%%
% plot KNN kinematic distributions standev image and relative colormap
close all;fNum = 1;
[~,currState,currKin] = plotKNNKinematicsSTD(f_orcoAll,meta.States2Plot_KNN,...
    meta.tSlice,'imagesc','relative',fNum);
plotKNNDistributions(currState,currKin,[plotFold subFoldName 'KNN_Relative_Heatmap/'],...
    '_relative_heatmap_allGenotypes_std')

%%
% plot turn optimality during inhibition for sharp turn
close all;fNum = 1;state = 1;
fNum = plotInhibitionOpt2(f_orcoAll,state,meta.tSlice2,fNum);
printFigures(fNum,[plotFold subFoldName 'Inhibition/'],'SharpTurn_TurnOptimality_inhibition_allGenotypes')

%%
% plot turn optimality during inhibition for curved walks
close all;fNum = 1;state = 2;
fNum = plotInhibitionOpt2(f_orcoAll,state,meta.tSlice2,fNum);
printFigures(fNum,[plotFold subFoldName 'Inhibition/'],'CurvedWalk_TurnOptimality_inhibition_allGenotypes')

%%
% plot turn optimality imagesc and absolute colormap
close all;fNum = 1;
[fNum,currState] = plotTurnOpt2(f_orcoAll,meta.States2Plot_Opt,meta.tSlice,...
    'imagesc','absolute',true,false,fNum);
currKin = cellstr(repmat('TurnOptimality',numel(currState),1))';
plotKNNDistributions(currState,currKin,[plotFold subFoldName 'KNN_Absolute_Heatmap/'],'_absolute_heatmap')

%%
% plot turn optimality imagesc and relative colormap
close all;fNum = 1;
[fNum,currState] = plotTurnOpt2(f_orcoAll,meta.States2Plot_Opt,meta.tSlice,...
    'imagesc','relative',true,false,fNum);
currKin = cellstr(repmat('TurnOptimality',numel(currState),1))';
h2 = findall(groot,'Type','figure');
plotKNNDistributions(currState,currKin,[plotFold subFoldName 'KNN_Relative_Heatmap/'],'_relative_heatmap')
printFigures(setdiff([h2.Number],1:numel(currKin)),[plotFold subFoldName],'KNN_TurnOptimality_OLS_dF')

%%
close all
end

function [] = plotKNNDistributions(currState,currKin,folderName,figExtension)
uState = unique(currState);
uKin = unique(currKin);
for i = 1:numel(uState)
    for j = 1:numel(uKin)
        figNum = find(strcmp(currState, uState(i)) & strcmp(currKin, uKin(j)));
        if ~isempty(figNum)
            figNum = [figNum nan];
            printFigures(figNum,folderName,...
                ['KNN_' uState{i} '_' uKin{j} figExtension])
        end
    end
end
end


function [] = printFigures(fNum,folder,figTitle)
if ~exist(folder, 'dir')
    mkdir(folder)
end

if numel(fNum)==1
    fRange = 1:max(fNum-1,1);
else
    fRange = fNum;
end
fRange(isnan(fRange)) = [];
if exist([folder figTitle '.ps'], 'file')==2
  delete([folder figTitle '.ps']);
end

for f = fRange
    figure(f);
    print('-bestfit','-painters','-dpsc2',[folder figTitle '.ps'],'-loose','-append');
end
ps2pdf('psfile', [folder figTitle '.ps'], 'pdffile', [folder figTitle '.pdf'], 'gspapersize', 'letter',...
    'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
    'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
    'gslibpath','C:\Program Files\gs\gs9.50\lib');
end

%%% some old code
% % plot spatial-temporal position
% close all;fNum = 1;
% for i = 1:numel(f_orcoAll)
%     if mod(i,6)==1
%         figure(fNum);set(gcf,'Position',[2 42 1000 924])
%         fNum = fNum+1;k = 1;
%     end
%     subplot(3,2,k);
%     [y_emp,xCent,yCent] = plotRadialPosition(f_orcoAll{i},meta.border,'after',true,...
%         'dt',2.*f_orcoAll{i}.fs,'clims',[0 0.01]);toc;
%     title(f_orcoAll{i}.id)
%     k = k+1;
% end
% printFigures(fNum,'Figures/','SpatialTemporal_Density_AllGenotypes')
%%%