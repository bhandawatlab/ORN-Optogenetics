function [] = EnsureReqFiles(meta)

if ~exist(meta.DestPath, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.DestPath]);
else
    disp([meta.DestPath " exists"]);
end

if ~exist(meta.foldDataModel, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.foldDataModel ]); 
else
    disp([meta.foldDataModel " exists"]);  
end

if ~exist(meta.plotFold, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.plotFold]);
else
    disp([meta.plotFold " exists"]);
end

if ~exist(meta.foldRaw, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.foldRaw]);
else
    disp([meta.foldRaw " exists"]);
end

if ~exist(meta.foldStim, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.foldStim]);
else
    disp([meta.foldStim " exists"]);
end

if ~exist(meta.folderData, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.folderData]);
else
    disp([meta.folderData " exists"]);
end

if ~exist(meta.foldSpk, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.foldSpk]);
else
    disp([meta.foldSpk " exists"]);
end

if ~exist(meta.LFPFilterFile, 'file')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.LFPFilterFile]);
else
    disp([meta.LFPFilterFile " exists"]);
end

if ~exist(meta.RateFilterFile, 'file')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.RateFilterFile]);
else
    disp([meta.RateFilterFile " exists"]);
end

if ~exist(meta.spikingDataFile, 'file')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.spikingDataFile]);
else
    disp([meta.spikingDataFile " exists"]);
end

if ~exist(meta.trainingDataFolder, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.trainingDataFolder]);
else
    disp([meta.trainingDataFolder " exists"]);
end

if ~exist(meta.syntheticFlyFold, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.syntheticFlyFold]);
else
    disp([meta.syntheticFlyFold " exists"]);
end

if ~exist(meta.syntheticPlotFold, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.syntheticPlotFold]);
else
    disp([meta.syntheticPlotFold " exists"]);
end

if ~exist(meta.summationModelFold, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.summationModelFold]);
else
    disp([meta.summationModelFold " exists"]);
end

if ~exist(meta.calibrationFolder, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.calibrationFolder]);
else
    disp([meta.calibrationFolder " exists"]);
end

if ~exist(meta.calibrationFile, 'file')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.calibrationFile]);
else
    disp([meta.calibrationFile " exists"]);
end

if ~exist(meta.foldName, 'dir')
    disp(["Entered Path Doesn't Exist, please enter valid path in Config file" meta.foldName]);
else
    disp([meta.foldName " exists"]);
end

end