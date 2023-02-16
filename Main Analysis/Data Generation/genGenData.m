function [] = genGenData(genAll,meta)
% This function creates a general file for each genotype with information
% about all flies for that genotype
%
% Inputs:
%    genAll: cell array of genotypes to consider
%    meta: Type of points we want to return (3 options)
%       meta.foldRaw = folder with the original data files
%       meta.genFold = folder to save the general files to

% get directory names
dataFold = meta.foldRaw;
genFolder = meta.folderData;
warning('off','MATLAB:load:variableNotFound');
nGen = numel(genAll);

% get the names of each file in each genotype directory
currfiles = cell(1,nGen);
genotype = struct('name',[],'files',[]);
for i = 1:nGen
    gen = genAll{i};
    currfiles{i} = dir([dataFold '\' gen '\*.mat']);
    genotype.name{i} = gen;
    for j = 1:length(currfiles{i})
        genotype.files{i,j} = currfiles{i}(j).name;
    end
end

% main loop to process data
fly = cell(1,nGen);curvError = cell(1,nGen);xyError = cell(1,nGen);
for i = 1 :nGen
    nFly = length(currfiles{i});
    sAll = cell(1,nFly);
    sArenaAll = cell(1,nFly);
    for j = 1:nFly
        fName = genotype.files{i,j};
        oldFolder = [dataFold '\' genotype.name{i}];
        clearvars s sArena
        load([oldFolder '\' fName],'s','fs','sArena');
        if ~exist('s')
            clearvars s fs sArena
            load([oldFolder '\' fName],'data');
            
            sArena.rad = mean([range(data.xpts)./2,range(data.ypts)./2]);
            sArena.arenaCent(1) = (max(data.xpts)+min(data.xpts))./2;
            sArena.arenaCent(2) = (max(data.ypts)+min(data.ypts))./2;
            sArena.cF = 4./sArena.rad;
            
            s.Center.x = data.xpts';
            s.Center.y = data.ypts';
            s.Center.xUS = data.xpts';
            s.Center.yUS = data.ypts';
            s.Head.x = data.xpts';
            s.Head.y = data.ypts';
            
            s.AngVec = myatan(gradient(s.Head.x),gradient(s.Head.y),'degrees',2);
            s.Kinematics.thrust = data.speed1'.*10;
            s.Kinematics.slip = zeros(size(s.Kinematics.thrust));
            s.Kinematics.yaw = gradient(s.AngVec);
            
            s.LightOn = 5400;
            fs = meta.fs;
            save([oldFolder '\' fName],'data','s','fs','sArena');
        end
        sAll{j} = s; sArenaAll{j} = sArena;
        
        try
            % calculate the curvature
            [curv,xErr,yErr,err] = calcCurvature(s,meta.fs,sArena);            
        catch
            [curv,~] = CalcCurvature(s);
            xErr = nan;yErr = nan;err = nan;
        end
        % assign the values to cells
        curvError{i}(j) = err;
        xyError{i}(:,j) = [xErr;yErr];
        fly{i}.curvature{1,j} = curv;
    end
    
    % create file with all the track information for flies within the 
    % current genotype
    sAll = createDataGen(sAll,sArenaAll,fly,genotype,curvError,xyError,meta.fs,i,genFolder,meta);
%     for KK = 1:nFly
%         s = sAll{KK};
%         sArena = sArenaAll{KK};
%         fName = genotype.files{K,KK};
%         newFolder = [dataNewFolder '\' genotype.name{K}];
%         if ~exist(newFolder, 'dir')
%             mkdir(newFolder)
%         end
%         save([newFolder '\' fName],'s','fs','sArena');
%     end
end
end

function [curv,xErr,yErr,err] = calcCurvature(s,fs,sArena)
        
% check if we can create x,y,ang from the thrust, slip, yaw values.
[xErr,yErr] = TSYtoXY(s,fs,sArena);

% create curvature from xy positions
[curv,~] = CalcCurvature(s);

% check if we can create the true curvature using thrust, slip, yaw
[~,err,~] = CalcCurvatureTSY(s,fs,sArena);

end

function [sAll] = createDataGen(sAll,sArenaAll,fly,genotype,curvError,xyError,fs,K,genFolder,meta)
% This just an assignment block
for KK = 1:length(sAll)
    s = sAll{KK};sArena = sArenaAll{KK};
        
    if max(sqrt(s.Center.x.^2+s.Center.y.^2))>4
        [s] = convertToUnitCircle(s,sArena);
        s.Center.x = s.Center.x.*4;% hard coded 
        s.Center.y = s.Center.y.*4;
        s.Head.x = s.Head.x.*4;
        s.Head.y = s.Head.y.*4;
    else
        a = 1;
    end
    
    
    Data.xHead(KK,1:length(s.Head.x)) = s.Head.x;
    Data.yHead(KK,1:length(s.Head.y)) = s.Head.y;
    Data.x(KK,1:length(s.Center.x)) = s.Center.x;
    Data.y(KK,1:length(s.Center.y)) = s.Center.y;
    Data.thrust(KK,1:length(s.Kinematics.thrust)) = s.Kinematics.thrust;
    Data.slip(KK,1:length(s.Kinematics.slip)) = s.Kinematics.slip;
    Data.yaw(KK,1:length(s.Kinematics.yaw)) = s.Kinematics.yaw;
    Data.ang(KK,1:length(s.AngVec)) = s.AngVec;
    Data.curv(KK,1:length(fly{K}.curvature{KK})) = fly{K}.curvature{KK};
    
    Data.lightOn(KK) = s.LightOn;

    Arena.center(:,KK) = sArena.arenaCent';
    Arena.rad(KK) = sArena.rad;
    Arena.cF(KK) = sArena.cF;

    gen.files = genotype.files(K,:);
    gen.files = gen.files(~cellfun('isempty',gen.files));
    
    sAll{KK} = s;
    
    len_data(KK) = length(s.Head.x);
end

for KK = 1:length(sAll)
    Data.xHead(KK,len_data(KK)+1:end) = nan;
    Data.yHead(KK,len_data(KK)+1:end) = nan;
    Data.x(KK,len_data(KK)+1:end) = nan;
    Data.y(KK,len_data(KK)+1:end) = nan;
    Data.thrust(KK,len_data(KK):end) = nan;
    Data.slip(KK,len_data(KK):end) = nan;
    Data.yaw(KK,len_data(KK):end) = nan;
    Data.ang(KK,len_data(KK)+1:end) = nan;
    Data.curv(KK,len_data(KK):end) = nan;
end

gen.name = genotype.name{K};
Error.curv = curvError{K};
Error.xy = xyError{K};

% save the data
fileName  = [genFolder '\' gen.name '_' meta.d '.mat'];
if isfile(fileName)
     save(fileName,'Data','Arena','fs','Error','gen','-append','-v7.3');
else
     save(fileName,'Data','Arena','fs','Error','gen','-v7.3');
end

%save([genFolder '\' genotype.name{K} '_Nov13.mat'],'Data','Arena','fs','Error','gen','-v7.3');

end







