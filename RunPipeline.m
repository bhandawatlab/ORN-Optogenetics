clear all;
close all hidden
prompt = "Using Existing Experimental Data, Y or N:";
txt = input(prompt,"s");

%Initialization
if txt =="Y"
    prompt = "Enter Destination File Path:";
    Destination_path = input(prompt,"s");
    while ~exist(Destination_path, 'dir')
        Destination_path = input("Entered Destination Path Doesn't Exist, please enter valid path again: ","s");
    end
    [meta] = Initialize_Params(Destination_path);
    %EnsureReqFiles(meta)
else
    disp("Generating New Experimental Analysis Dataset Folder");
    dt = datestr(now, 'yyyymmdd-HH_MM');% Download Datset
    if ~exist('Dataset.zip', 'file')
        outfilename = websave('Dataset.zip','https://www.dropbox.com/s/qjyx5voz82onfic/Data.zip?dl=1');
    end
    EXP_folder = ['Experiment_' dt];
    if ~exist(EXP_folder, 'dir')
            mkdir(EXP_folder);
    end
    unzip('Dataset.zip',['Experiment_' dt]);
    disp("Dataset Downloaded and Extracted");

    Destination_path = [pwd '\' EXP_folder];
    [meta] = Initialize_Params(Destination_path);
end

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

%% generate consolidated data file from individual tracked data files
if meta.dataFromRaw
    %--------------------------------------------------------------------------
    % run this section when you are using data files from the tracking code
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % generate Data Files (new Data) and Data Gen (consolidated Data by genotype)
    disp('Generating new data files...');
    genGenData(genAll,meta);
    %--------------------------------------------------------------------------
%     % findbest fit parameters to delineate between sharp turn and curved walks
%     disp('Finding best parameters...')
%     [GlobMinX,GlobMinCFit] = getBestFit(meta);
    %--------------------------------------------------------------------------
    % generate sharp turn/curved walk data
    disp('Generating sharp turn curved walk data files...');
    load([Destination_path '/Data/BestFit5.mat'],'GlobMinX');
    genSTCWDat(genAll,meta,GlobMinX);
    %--------------------------------------------------------------------------
end

EnsureReqFiles(meta)
%% generate light intensity from radial positions
label = 'Head';
GenerateLightIntensity(label,meta);

%% get the spiking data
GenerateSpikingData(genAll,meta);

%% generate the empirical fly object
GenerateEmpiricalFlies(genAll,meta);

%% plot all figures
PlotAllFiguresEmp(genAll,gen_Retinal,meta)

%% ORN to speed/curv linear filter analysis
gen = gen_Retinal{1};%'Orco Retinal'
load(strcat(string(meta.foldDataModel),'\',gen,'_',meta.d,meta.ext,'.mat'),'f_orco');%_allTime
transitionOnly = true;
plotFigure = true;
tic;[~,~] = linearFilterAnalysisCW2TurnTransition(f_orco,transitionOnly,meta,plotFigure);toc;
transitionOnly = false;
tic;[~,~] = linearFilterAnalysisCW2TurnTransition(f_orco,transitionOnly,meta,plotFigure);toc;
% kinematics linear filters when entering and leaving the light border
[~,~,~,~] = linearFilterAnalysisSpeedCurvature(f_orco,meta,plotFigure);
linearFilterAnalysisSpeedCurvatureAllTime(f_orco,meta);
linearFilterStateKinematics(gen_Retinal{1},meta);
close all

%% generate synthetic flies
agentModelMeta.delay = 0; %time since peak of filter
agentModelMeta.dur = 0; %duration of crossing filters%2
if meta.adaptation == true
    % generate synthetic flies
    agentModelHandle(gen_Retinal,meta,agentModelMeta,true)
    % get the time averaged KNN space of these synthetic flies
    meta2 = meta;
    meta2.adaptation = false;
    meta2.tSlice = 0;
    meta2.tSlice2 = 0;
    meta2.zGrid = [0 180].*meta.fs;
    meta2.ratio = [30, 10, 1000.*meta.fs];
    meta2.ext = '_allTime';
    meta2.foldName = 'TimeAveraged/';
    agentModelHandle(gen_Retinal,meta2,agentModelMeta,false)
    PlotAllFiguresSynth(gen_Retinal,meta2,agentModelMeta.delay,agentModelMeta.dur(1));
else
    agentModelHandle(gen_Retinal,meta,agentModelMeta,true)
end

%% plot empirical synthetic comparisons
PlotAllFiguresSynth(gen_Retinal,meta,agentModelMeta.delay,agentModelMeta.dur(1));
PlotAllFiguresEmpSynth(gen_Retinal,meta,agentModelMeta);

%% plot rules of summation analysis
% gen_Retinal = {'Orco Retinal','Ir8a Retinal','Orco Ir8a Retinal', ...
%     'Or42a Retinal','Or42b Retinal', 'Or92a Retinal', 'Ir64a Retinal', 'Ir75a Retinal', ...
%     'Or42b Or92a Retinal', 'Ir64a Ir75a Retinal', 'Ir64a Or42b Retinal', ...
%     'Ir64a Ir75a Or42b Retinal', 'Or42a Or42b Or92a Retinal','Ir64a Ir75a Or42b Retinal'};
% gen_Control = {'Orco Control','Ir8a Control','Orco Ir8a Control', ...
%     'Or42a Control','Or42b Control', 'Or92a Control', 'Ir64a Control', 'Ir75a Control', ...
%     'Or42b Or92a Control', 'Ir64a Ir75a Control', 'Ir64a Or42b Control', ...
%     'Ir64a Ir75a Or42b Control', 'Or42a Or42b Or92a Control','Ir64a Ir75a Or42b Control'};
rulesOfSummationAnalysis(gen_Retinal,meta)

genTimeInside = {'Or42b Or92a Control','Or42b Or92a Retinal',...
    'Ir64a Ir75a Or42b Control','Ir64a Ir75a Or42b Retinal',};
plotTimeInsideLightZoneVsBorderRegion(genTimeInside,meta);
close all
