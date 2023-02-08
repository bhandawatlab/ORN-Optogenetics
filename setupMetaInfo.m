function [meta] = setupMetaInfo(adaptation,plotFig,plotSupplements)
% setupMetaInfo  Function to set up meta information. This includes:
%       1. add relevant path
%       2. setup folders for data, analysis, and plotting
%       3. parameters regarding movement states, arena, KNN analysis, etc
%       
%   Inputs: adaptation = true/false for time adaptation
%           plotFig = true/false to plot figure output
%           plotSupplements = true/false to plot supplemental analysis
%           (e.g. those that are not genotype specific)
%   

% add paths
addpath(genpath([pwd '/Util']))
addpath(genpath([pwd '/Subroutines']))
addpath(genpath([pwd '/Main Analysis']))
addpath(genpath([pwd '/PlottingFunctions']))
addpath(genpath([pwd '/Model']))

% date identifier
meta.d = 'March2022';

% set folders and supporting files
meta.foldRaw = [pwd '\Data\DataRaw'];
meta.folderData = [pwd '\Data\DataGen'];
meta.foldStim = [pwd '\Data\DataStim'];
meta.foldSpk = [pwd '\Data\DataSpike'];
meta.folderObject = [pwd '\DataModel'];
meta.LFPFilterFile = [pwd '\Data\LFP_2_25.mat'];
meta.RateFilterFile = [pwd '\Data\LFP2Rate_4_12.mat'];
meta.Intensity2VoltageFile = [pwd '\Data\Intensity_space2.mat'];
meta.spikingDataFile = [pwd '\Data\empdataNew.mat'];
meta.trainingDataFolder = [pwd '\Data\Training Data\'];
meta.syntheticFlyFold = [pwd '\DataRT\'];
meta.summationModelFold = [pwd '\Data\SummationModel\'];
meta.calibrationFolder = [pwd '\Data\Calibration'];
meta.calibrationFile = '\LIR_exp_1.xlsx';

meta.plotFold = [pwd '\Figures\'];
meta.syntheticPlotFold = [meta.plotFold '\EmpSynth\'];

% create any folder that does not exist
createFolders(meta);

meta.fs_linFilter = 100;% linear filter frame rate

% generating flies parameters
meta.border = 1.25;% cm
meta.fs = 30;% frames/s
meta.stopThresh = 0.5;% mm/s
meta.rBound = 4;% cm
meta.saveData = true;
meta.adaptation = adaptation;

% plotting parameters
meta.plotFig = plotFig;%
meta.plotSupplements = plotSupplements;%false

%% KNN parameters
% time period to average over before start of each state trajectory
meta.timeInterval = [0 200];% in ms
% time since first entry slices shown in the paper
if adaptation == true
    meta.tSlice = [0:15:115];meta.tSlice2 = [0:30:180]; 
else
    meta.tSlice = 0;meta.tSlice2 = 0;
end
meta.States2Plot_Inh = 1:4;% sharp turn and curved walk speed/curvature
meta.States2Plot_KNN = 1:8;% sharp turn, curved walk, and stop parameters
meta.States2Plot_Opt = 1:3;% ST, CW, and stop turn optimality

% set the time slice and radii of the bounding ellipse for KNN estimations
if adaptation == true
    %----------------------------------------------------------------------
    % in slices of 5 seconds time since first entry
    meta.zGrid = [0:5:180].*meta.fs;
    meta.ratio = [30, 10, 20.*meta.fs];
    meta.ext = '';
    meta.foldName = 'Adaptation/';
    %----------------------------------------------------------------------
else
    % no time since first entry
    meta.zGrid = [0 180].*meta.fs;
    meta.ratio = [30, 10, 1000.*meta.fs];
    meta.ext = '_allTime';
    meta.foldName = 'TimeAveraged/';
    %----------------------------------------------------------------------
end
meta.xGrid = -150:15:150;
meta.yGrid = 0:1:55;
meta.Thresh = [];
meta.K = [];

end


function [] = createFolders(meta)
% set up what folders to create
fold{1} = meta.foldRaw;
fold{2} = meta.folderData;
fold{3} = meta.foldStim;
fold{4} = meta.foldSpk;
fold{5} = meta.folderObject;
fold{6} = meta.trainingDataFolder;
fold{7} = meta.syntheticFlyFold;
fold{8} = meta.syntheticPlotFold;
fold{9} = meta.summationModelFold;

% create the necessary folder to house the data
for i = 1:numel(fold)
    if ~exist(fold{i}, 'dir')
        mkdir(fold{i})
    end
end
end












