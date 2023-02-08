clear;close all hidden
%%
adaptation = false;
plotFig = false;
plotSupplements = true;
[meta] = setupMetaInfo(adaptation,plotFig,plotSupplements);

% set genotypes
gen_Retinal = {'Orco Retinal','Ir8a Retinal','Orco Ir8a Retinal',...
    'Or42a Retinal','Or42b Retinal', 'Or92a Retinal', 'Ir64a Retinal', 'Ir75a Retinal',...
    'Or42b Or92a Retinal', 'Ir64a Ir75a Retinal', 'Ir64a Or42b Retinal', ...
    'Ir64a Ir75a Or42b Retinal', 'Or42a Or42b Or92a Retinal','Ir64a Ir75a Or42b New Retinal'};
gen_Control = {'Orco Control','Ir8a Control','Orco Ir8a Control',...
    'Or42a Control','Or42b Control', 'Or92a Control', 'Ir64a Control', 'Ir75a Control',...
    'Or42b Or92a Control', 'Ir64a Ir75a Control', 'Ir64a Or42b Control', ...
    'Ir64a Ir75a Or42b Control', 'Or42a Or42b Or92a Control','Ir64a Ir75a Or42b Control'};
genAll = reshape([gen_Retinal;gen_Control],[],1);

%% generate consolidated data file from individual tracked data files
if meta.dataFromRaw
    %--------------------------------------------------------------------------
    % run this section when you are using data files from the tracking code
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % generate Data Files (new Data) and Data Gen (consolidated Data by genotype)
    disp('Generating new data files...')
    genGenData(genAll,meta);
    %--------------------------------------------------------------------------
    % findbest fit parameters to delineate between sharp turn and curved walks
    disp('Finding best parameters...')
    [GlobMinX,GlobMinCFit] = getBestFit(meta);
    %--------------------------------------------------------------------------
    % generate sharp turn/curved walk data
    disp('Generating sharp turn curved walk data files...')
    load('Data/BestFit5.mat','GlobMinX');
    genSTCWDat(genAll,meta,GlobMinX);
    %--------------------------------------------------------------------------
end
%% generate light intensity from radial positions
label = 'Head';
GenerateLightIntensity(label,meta);

%% get the spiking data
GenerateSpikingData(meta);

%% generate the empirical fly object
GenerateEmpiricalFlies(genAll,meta);

%% plot all figures
PlotAllFiguresEmp(genAll,gen_Retinal,meta)

%% ORN to speed/curv linear filter analysis
gen = gen_Retinal{1};%'Orco Retinal'
load([meta.folderObject '\' gen '_' meta.d meta.ext '.mat'],'f_orco');%_allTime
transitionOnly = true;
plotFigure = true;
[~,~] = linearFilterAnalysisCW2TurnTransition(f_orco,transitionOnly,meta,plotFigure);
transitionOnly = false;
[~,~] = linearFilterAnalysisCW2TurnTransition(f_orco,transitionOnly,meta,plotFigure);
% kinematics linear filters when entering and leaving the light border
[~,~,~,~] = linearFilterAnalysisSpeedCurvature(f_orco,meta,plotFigure);
linearFilterAnalysisSpeedCurvatureAllTime(f_orco,meta);

%% plot rules of summation analysis
rulesOfSummationAnalysis(gen_Retinal,meta)

%% generate synthetic flies
agentModelMeta.delay = 0;% time since peak of filter
agentModelMeta.dur = [0];% duration of crossing filters%2
if meta.adaptation == true
    agentModelHandle(gen_Retinal,meta,agentModelMeta,true)
else
    agentModelHandle(gen_Retinal,meta,agentModelMeta,true)
end

%% plot empirical synthetic comparisons
%PlotAllFiguresSynth(gen_Retinal,meta,agentModelMeta.delay,agentModelMeta.dur);
PlotAllFiguresEmpSynth(gen_Retinal,meta,agentModelMeta);

%% plot rules of summation analysis
rulesOfSummationAnalysis(gen_Retinal,meta)

%% svd analysis
load([meta.folderObject '\' gen '_' meta.d meta.ext '.mat'],'f_orco');%_allTime
figureFile = 'StateSpaceEmbeddingDimensionsFNN_200msDelay';
SVD_embeddingAnalysis(f_orco,meta,figureFile);

%% ORN to speed/curv linear filter analysis
gen = gen_Retinal{1};%'Orco Retinal'
load([meta.folderObject '\' gen '_' meta.d meta.ext '.mat'],'f_orco');%_allTime
transitionOnly = true;
plotFigure = true;
[~,~] = linearFilterAnalysisCW2TurnTransition(f_orco,transitionOnly,meta,plotFigure);
transitionOnly = false;
[~,~] = linearFilterAnalysisCW2TurnTransition(f_orco,transitionOnly,meta,plotFigure);
% kinematics linear filters when entering and leaving the light border
[~,~,~,~] = linearFilterAnalysisSpeedCurvature(f_orco,meta,plotFigure);
linearFilterAnalysisSpeedCurvatureAllTime(f_orco,meta);


















