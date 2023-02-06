function [GlobMinX,GlobMinCFit] = getBestFit(genAll,meta)
% This function fits 6 parameters to delineate between sharp turns and 
% curved walks
%
% Inputs:
%    opts:
%       opts.folderData = folder with fly tracks
% 
% Outputs:
%    GlobMinX: The best fit parameters
%    GlobMinCFit: RMSE of the fit

% get data folder
genFolder = meta.folderData;
nGen = numel(genAll);

Ndx = 1;s = {};curv = {};
for i = 1:nGen
    fileName  = [genFolder '\' genAll{i} '_' meta.d '.mat'];
    load(fileName,'Data');
    nFly = numel(Data.lightOn);
    
    % loop through each fly for the genotype
    for j = 1:nFly
        curvTmp = Data.curv(j,:);
        
        % calulate the speed and radial position
        spd = sqrt(Data.thrust(j,:).^2+Data.slip(j,:).^2);
        r = sqrt(Data.x(j,:).^2+Data.y(j,:).^2);
        % assign stops and boundary conditions
        stop = spd<meta.stopThresh;
        boundary = r(:,1:end-1)>meta.rBound;
        stop(boundary) = false;
        
        % assign short (<3 frame) tracks not assigned to stop/boundary to 
        % stop/boundary
        [startNdx,endNdx,type] = startEndSeq(stop|boundary);
        startNdx = startNdx(type == 0);
        endNdx = endNdx(type == 0);
        shortTracks = find(endNdx-startNdx+1<3);
        if sum(shortTracks>0)
            for k = 1:length(shortTracks)
                if stop(max(startNdx(shortTracks(k))-1,1))
                    stop(startNdx(shortTracks(k)):endNdx(shortTracks(k)))=true;
                else
                    boundary(startNdx(shortTracks(k)):endNdx(shortTracks(k)))=true;
                end
            end
        end
        startNdx(shortTracks) = [];
        endNdx(shortTracks) = [];
        
        % assign cell arrays of x, y, and curvature trajectories to use in
        % the fitting process
        for k = 1:length(startNdx)
            s{Ndx}.x = Data.x(startNdx(k):endNdx(k));
            s{Ndx}.y = Data.y(startNdx(k):endNdx(k));
            curv{Ndx} = curvTmp(startNdx(k):endNdx(k)-1);
            Ndx = Ndx+1;
        end
        
    end
end
% fit to the trajectories
[GlobMinX,GlobMinCFit] = calcFit(s,curv);
% save the best fit parameters
save('BestFit5','GlobMinX','GlobMinCFit','-v7.3')
end

function [GlobMinX,GlobMinCFit] = calcFit(s,curv)
% do not display plots
pltInfo = [];
% set up objective function
obj = @(x)objFunc(x(1),x(2),x(3),x(4),x(5),x(6),s,curv,pltInfo);
gs = GlobalSearch('Display','iter','NumTrialPoints',2000);

% fit the objective function
tic;x = objFuncWrap(obj,gs);toc

% refit the objective function until either 5 iterations has passed or a
% "low" error has been acheived
j = 0;
while (obj(x)>5 && j<5)
    display(j)
    xTmp = objFuncWrap(obj,gs);
    % keep the new parameters if it fit better than the previous best set
    % of parameters
    if obj(xTmp) < obj(x)
        x = xTmp;
    end
    j = j+1;
end
toc
GlobMinX = x;
% get the fit from the "best" parameter set
[GlobMinCFit,~,~] = obj(x);
end


