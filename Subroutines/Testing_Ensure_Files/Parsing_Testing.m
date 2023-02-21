function [meta] = Initialize_Params()    % clc
    % clear all;
    % close all;
    % Dataset = 'Dataset'
    % if ~exist(Dataset, 'dir')
    %         mkdir(Dataset);
    % end


    %fin = fopen('Test.xlsx');
    %disp(data)
    %A = csvread('Test.csv') %, 'columns', {'X, 'Y', 'Z'});
    %disp(A)
    A = readcell(['Test.csv']); %textscan(fin,'%s','Delimiter','\n')
    B = string(A(:, 1));

    dateValIdx = find(B=='d')
    meta.d = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='foldDataModel')
    meta.foldDataModel = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='foldRaw')
    meta.foldRaw = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='foldStim')
    meta.foldStim = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='folderData')
    meta.folderData = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='foldSpk')
    meta.foldSpk = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='LFPFilterFile')
    meta.LFPFilterFile = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='RateFilterFile')
    meta.RateFilterFile = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='Intensity2VoltageFile')
    meta.Intensity2VoltageFile = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='spikingDataFile')
    meta.spikingDataFile = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='trainingDataFolder')
    meta.trainingDataFolder = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='syntheticFlyFold')
    meta.syntheticFlyFold = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='syntheticPlotFold')
    meta.syntheticPlotFold = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='summationModelFold')
    meta.summationModelFold = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='calibrationFolder')
    meta.calibrationFolder = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='calibrationFile')
    meta.calibrationFile = [strcat(pwd,erase(A(dateValIdx,2),"'"))];

    dateValIdx = find(B=='fs_linFilter')
    meta.fs_linFilter = str2num(string(A(dateValIdx,2)));

    dateValIdx = find(B=='border')
    meta.border = str2num(string(A(dateValIdx,2)));

    dateValIdx = find(B=='fs')
    meta.fs = str2num(string(A(dateValIdx,2)));

    dateValIdx = find(B=='stopThresh')
    meta.stopThresh = str2num(string(A(dateValIdx,2)));

    dateValIdx = find(B=='rBound')
    meta.rBound = str2num(string(A(dateValIdx,2)));


    % create any folder that does not exist
    createFolders(meta);

    dateValIdx = find(B=='saveData')
    if str2num(string(A(dateValIdx,2)))==1
        meta.saveData = true
    else
        meta.saveData = false
    end
    dateValIdx = find(B=='adaptation')
    if str2num(string(A(dateValIdx,2)))==1
        meta.adaptation = true
    else
        meta.adaptation = false
    end

    dateValIdx = find(B=='plotFig')
    if str2num(string(A(dateValIdx,2)))==1
        meta.plotFig = true
    else
        meta.plotFig = false
    end

    dateValIdx = find(B=='plotSupplements')
    if str2num(string(A(dateValIdx,2)))==1
        meta.plotSupplements = true
    else
        meta.plotSupplements = false
    end
    dateValIdx = find(B=='timeInterval')
    meta.timeInterval = [str2num(string(A(dateValIdx,2))) str2num(string(A(dateValIdx,3)))];% in ms

    dateValIdx = find(B=='States2Plot_Inh')
    meta.States2Plot_Inh = str2num(string(A(dateValIdx,2))):str2num(string(A(dateValIdx,3)))
    dateValIdx = find(B=='States2Plot_KNN')
    meta.States2Plot_KNN = str2num(string(A(dateValIdx,2))):str2num(string(A(dateValIdx,3)))

    dateValIdx = find(B=='States2Plot_Opt')
    meta.States2Plot_Opt = str2num(string(A(dateValIdx,2))):str2num(string(A(dateValIdx,3)))


    if meta.adaptation == true
        meta.tSlice = [0:15:115];
        meta.tSlice2 = [0:30:180]; 
    else
        meta.tSlice = 0;
        meta.tSlice2 = 0;
    end

    if meta.adaptation == true
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

    dateValIdx = find(B=='xGrid')
    meta.xGrid = str2num(string(A(dateValIdx,2))):str2num(string(A(dateValIdx,3))):str2num(string(A(dateValIdx,4)));

    dateValIdx = find(B=='yGrid')
    meta.yGrid = str2num(string(A(dateValIdx,2))):str2num(string(A(dateValIdx,3))):str2num(string(A(dateValIdx,4)));

    meta.Thresh = [];
    meta.K = [];

    function [] = createFolders(meta)
    % set up what folders to create
    fold(1) = meta.foldRaw;
    fold(2) = meta.foldStim;
    fold(3) = meta.folderData;
    fold(4) = meta.foldSpk;
    fold(5) = meta.trainingDataFolder;
    fold(6) = meta.syntheticFlyFold;
    fold(7) = meta.syntheticPlotFold;
    fold(8) = meta.summationModelFold;
    fold(9) = meta.foldDataModel;


    % create the necessary folder to house the data
    for i = 1:numel(fold)
        if ~exist(string(fold(i)), 'dir')
            mkdir(string(fold(i)))
        else
            disp("Folder Already Exists")
        end
    end