function PlotAllFiguresSynth(genAll,meta,currDelay,currDur)

cond = ['_CWdur' strrep(num2str(currDur),'.','') ...
                '_CWdelay' strrep(num2str(currDelay),'.','')];

plotFolder = [meta.syntheticPlotFold];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting general analysis figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (for retinal only)
f_orcoAll = cell(1,numel(genAll));
for i = 1:numel(genAll)
    gen = genAll{i};
    fName = ['RT_run' gen '_' meta.d meta.ext cond '_flies.mat'];
    try
        load([meta.syntheticFlyFold fName],'synth_orco')
        f_orcoAll{i} = synth_orco;
    catch
        f_orcoAll{i} = [];
    end
end

%%
% create subfolders for plots
KNN_folders = {'KNN_Absolute_Heatmap','KNN_Relative_Heatmap','Inhibition'};
foldName = meta.foldName;
for i = 1:numel(KNN_folders)
    if ~exist([plotFolder foldName KNN_folders{i}], 'dir')
        mkdir([plotFolder foldName KNN_folders{i}])
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting general analysis figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% plot transition probability image and relative colormap
close all;fNum = 1;
[fNum] = plotTransitionProbability2(f_orcoAll,'imagesc','relative',fNum);
printFigures(fNum,[plotFolder foldName],'TransitionProbability_relative_heatmap_allGenotypes')

%%
% plot transition probability image and absolute colormap
close all;fNum = 1;
[fNum] = plotTransitionProbability2(f_orcoAll,'imagesc','absolute',fNum);
printFigures(fNum,[plotFolder foldName],'TransitionProbability_absolute_heatmap_allGenotypes')

%%
% plot baseline kinematic comparisons
% close all;fNum = 1;
% [fNum] = plotBaselineKinematics(f_orcoAll,meta.border,fNum);
% printFigures(fNum,[plotFolder foldName],'BaselineKinematics')

%%
% plot spatial probability map
close all;fNum = 1;
[fNum] = plotROISpatialProb(f_orcoAll,[0.1 0.1],'H',meta.border,fNum);
printFigures(fNum,[plotFolder foldName],'SpatialProbMap')

%
% plot inhibition kinematics
close all;fNum = 1;
[~,currState,currKin] = plotInhibitionKinematics2(f_orcoAll,...
    meta.States2Plot_Inh,meta.tSlice,fNum);
plotKNNDistributions(currState,currKin,[plotFolder foldName 'Inhibition/'],...
    '_inhibition_allGenotypes')

%%
% plot KNN kinematic distributions image and absolute colormap
close all;fNum = 1;
[~,currState,currKin] = plotKNNKinematics2(f_orcoAll,meta.States2Plot_KNN,...
    meta.tSlice,'imagesc','absolute',fNum);
plotKNNDistributions(currState,currKin,[plotFolder foldName 'KNN_Absolute_Heatmap/'],...
    '_absolute_heatmap_allGenotypes')

%%
% plot KNN kinematic distributions image and relative colormap
close all;fNum = 1;
[~,currState,currKin] = plotKNNKinematics2(f_orcoAll,meta.States2Plot_KNN,...
    meta.tSlice,'imagesc','relative',fNum);
plotKNNDistributions(currState,currKin,[plotFolder foldName 'KNN_Relative_Heatmap/'],...
    '_relative_heatmap_allGenotypes')

%%
% plot turn optimality during inhibition for sharp turn
close all;fNum = 1;state = 1;
fNum = plotInhibitionOpt2(f_orcoAll,state,meta.tSlice2,fNum);
printFigures(fNum,[plotFolder foldName 'Inhibition/'],'SharpTurn_TurnOptimality_inhibition_allGenotypes')

%%
% plot turn optimality during inhibition for curved walks
close all;fNum = 1;state = 2;
fNum = plotInhibitionOpt2(f_orcoAll,state,meta.tSlice2,fNum);
printFigures(fNum,[plotFolder foldName 'Inhibition/'],'CurvedWalk_TurnOptimality_inhibition_allGenotypes')

%%
% plot turn optimality imagesc and absolute colormap
close all;fNum = 1;
[fNum,currState] = plotTurnOpt2(f_orcoAll,meta.States2Plot_Opt,meta.tSlice,...
    'imagesc','absolute',true,false,fNum);
currKin = cellstr(repmat('TurnOptimality',numel(currState),1))';
plotKNNDistributions(currState,currKin,[plotFolder foldName 'KNN_Absolute_Heatmap/'],'_absolute_heatmap')

%%
% plot turn optimality imagesc and relative colormap
close all;fNum = 1;
[fNum,currState] = plotTurnOpt2(f_orcoAll,meta.States2Plot_Opt,meta.tSlice,...
    'imagesc','relative',true,false,fNum);
currKin = cellstr(repmat('TurnOptimality',numel(currState),1))';
h2 = findall(groot,'Type','figure');
plotKNNDistributions(currState,currKin,[plotFolder foldName 'KNN_Relative_Heatmap/'],'_relative_heatmap')
printFigures(setdiff([h2.Number],1:numel(currKin)),[plotFolder foldName],'KNN_TurnOptimality_OLS_dF')

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
