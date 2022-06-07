clear;close all hidden
%%
% add paths
addpath(genpath([pwd '/Util']))
addpath(genpath([pwd '/Subroutines']))
addpath(genpath([pwd '/Main Analysis']))
addpath(genpath([pwd '/PlottingFunctions']))
addpath(genpath([pwd '/Model']))

% date identifier
meta.d = 'April2022';

% set folders and supporting files
meta.foldStim = [pwd '\Data\DataStim'];
meta.foldSpk = [pwd '\Data\DataSpike'];
meta.LFPFilterFile = [pwd '\Data\LFP_2_25.mat'];
meta.RateFilterFile = [pwd '\Data\LFP2Rate_4_12.mat'];
meta.Intensity2VoltageFile = [pwd '\Data\Intensity_space2.mat'];
meta.spikingDataFile = [pwd '\Data\empdataNew.mat'];
meta.trainingDataFolder = [pwd '\Data\Training Data\'];
meta.syntheticFlyFold = [pwd '\DataRT\'];
meta.syntheticPlotFold = [pwd '\Figures\EmpSynth\'];
meta.summationModelFold = [pwd '\Data\SummationModel\'];

% set whether to generate flies with adaptation
adaptation = true;%false

% set genotypes
gen_Retinal = {'Orco Retinal','Ir8a Retinal','Orco Ir8a Retinal',...
    'Or42a Retinal','Or42b Retinal', 'Or92a Retinal', 'Ir64a Retinal', 'Ir75a Retinal',...
    'Or42b Or92a Retinal', 'Ir64a Ir75a Retinal', 'Ir64a Or42b Retinal', ...
    'Ir64a Ir75a Or42b Retinal', 'Or42a Or42b Or92a Retinal'};

gen_Control = {'Orco Control','Ir8a Control','Orco Ir8a Control',...
    'Or42a Control','Or42b Control', 'Or92a Control', 'Ir64a Control', 'Ir75a Control',...
    'Or42b Or92a Control', 'Ir64a Ir75a Control', 'Ir64a Or42b Control', ...
    'Ir64a Ir75a Or42b Control', 'Or42a Or42b Or92a Control'};

genAll = reshape([gen_Retinal;gen_Control],[],1);

% generating flies parameters
meta.border = 1.25;% cm
meta.fs = 30;% frames/s
meta.stopThresh = 0.5;% mm/s
meta.rBound = 4;% cm
meta.plotFig = false;%

% plotting parameters
meta.plotSupplements = false;%false
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
    meta.foldName = 'Habituation/';
    %----------------------------------------------------------------------
else
    % no time since first entry
    meta.zGrid = [0 180].*meta.fs;
    meta.ratio = [30, 10, 1000.*meta.fs];
    meta.ext = '_allTime';
    meta.foldName = 'TimeAveraged/';
    %----------------------------------------------------------------------
end

%% get the spiking data
GenerateSpikingData(meta)

%% generate the empirical fly object
GenerateEmpiricalFlies(genAll,meta)

%% plot all figures
PlotAllFiguresEmp(genAll,gen_Retinal,meta)

%% generate synthetic flies
agentModelMeta.tau = 0.5;%[50]./100;%[0.1 0.5 1 2];0.1
agentModelMeta.C = 0;%0.5;%[0.1 0.5 1 2];0.1
if adaptation == true
    agentModelHandle(gen_Retinal,meta,agentModelMeta,true)
else
    agentModelHandle(gen_Retinal,meta,agentModelMeta,false)
end

%% plot empirical synthetic comparisons
PlotAllFiguresEmpSynth(gen_Retinal,meta,agentModelMeta)
PlotAllFiguresSynth(gen_Retinal,meta,agentModelMeta.tau,agentModelMeta.C);

%% plot rules of summation analysis
rulesOfSummationAnalysis(gen_Retinal,meta)




















