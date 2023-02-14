function [] = genSTCWDat(genAll,meta,GlobMin)
% This function delineate tracks between sharp turns and curved walks
%
% Inputs:
%    genAll: cell array of genotypes to consider
%    GlobMin: structure of parameters to delineate sharp turn and curved
%       walks
%    meta: Type of points we want to return (3 options)
%       meta.dataFold = folder with the data

% get data folder
genFolder = meta.folderData;
nGen = numel(genAll);

errAll = cell(1,nGen);
for i = 1:nGen
    fileName  = [genFolder '\' genAll{i} '_' meta.d '.mat'];
    load(fileName,'Data');
    nFly = numel(Data.lightOn);
    
    figNum = 1;%k = 1;
    curvPks = [];curvWalks = [];stopCond = [];boundCond = [];
    err = zeros(nFly,1);
    % loop through each fly
    for j = 1:nFly
        curv = Data.curv(j,:);
        
        % calulate the speed and radial position
        spd = sqrt(Data.thrust(j,:).^2+Data.slip(j,:).^2);
        r = sqrt(Data.x(j,:).^2+Data.y(j,:).^2);
        % assign stops and boundary conditions
        stop = spd<meta.stopThresh;
        boundary = r(:,1:end-1)>(meta.rBound-0.15);% boundary if fly body is 1.5 mm (0.15 cm) away from the arena boundary
        stop(boundary) = false;
        
        % assign short (<3 frame) tracks not assigned to stop/boundary to
        % stop/boundary
        [startNdx,endNdx,type] = startEndSeq(stop|boundary);
        startNdx = startNdx(type == 0);
        endNdx = endNdx(type == 0);
        shortTracks = find(endNdx-startNdx+1<3);
        if sum(shortTracks>0)>0
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
        
        % get the stop and boundary condition structures
        [stopCond] = computeStopBorder(stop,Data.x(j,:),Data.y(j,:),curv,j,stopCond);
        [boundCond] = computeStopBorder(boundary,Data.x(j,:),Data.y(j,:),curv,j,boundCond);
        
        % initiate cells
        curvPks.max{j} = []; curvPks.tot{j} = []; curvPks.dur{j} = []; curvPks.avg{j} = [];
        curvPks.ndx{j} = []; curvPks.all{j} = []; curvPks.dirRelCenterBef{j} = [];
        curvPks.dirRelCenterAft{j} = [];
        curvWalks.max{j} = []; curvWalks.tot{j} = []; curvWalks.dur{j} = []; curvWalks.avg{j} = [];
        curvWalks.ndx{j} = []; curvWalks.all{j} = []; curvWalks.dirRelCenterBef{j} = [];
        curvWalks.dirRelCenterAft{j} = [];
        
        % update the figure conditions
        if mod(j,10) == 1
            figure(figNum);set(gcf,'Position',[2 42 638 600])
            p = 1;figNum = figNum+1;
        end
        
        % plot the entire fly trajectory
        figure(figNum-1);subplot(5,4,p);
        plot3(Data.x(j,:),Data.y(j,:),1:1:length(Data.x(j,:)),'g','LineWidth',0.5);
        hold on;view(2)
        
        % loop through each non boundary/stop track
        for k = 1:length(startNdx)
            % initialize x, y, and curvature trajectories to use in the
            % fitting process
            s2{1}.x = Data.x(j,startNdx(k):endNdx(k)+1);
            s2{1}.y = Data.y(j,startNdx(k):endNdx(k)+1);
            curv2{1} = curv(startNdx(k):endNdx(k));
            % give the subplot parameters
            pltInf = [5,4,p];
            % set up objective function
            obj = @(x)objFunc(x(1),x(2),x(3),x(4),x(5),x(6),s2,curv2,pltInf);
            [err(j,k),curvPksTmp,curvWalksTmp] = obj(GlobMin);
            
            % assign fields for sharp turns and curved walks
            for KKKK = 1:length(curvPksTmp.all)
                curvPksTmp.all{KKKK}(:,2) = curvPksTmp.all{KKKK}(:,2)+startNdx(k)-1;
            end
            for KKKK = 1:length(curvWalksTmp.all)
                curvWalksTmp.all{KKKK}(:,2) = curvWalksTmp.all{KKKK}(:,2)+startNdx(k)-1;
            end
            curvPks.max{j} = [curvPks.max{j} curvPksTmp.max];
            curvPks.tot{j} = [curvPks.tot{j} curvPksTmp.tot];
            curvPks.dur{j} = [curvPks.dur{j} curvPksTmp.dur];
            curvPks.avg{j} = [curvPks.avg{j} curvPksTmp.avg];
            curvPks.ndx{j} = [curvPks.ndx{j} curvPksTmp.ndx+startNdx(k)-1];
            curvPks.all{j} = [curvPks.all{j} curvPksTmp.all];
            curvPks.dirRelCenterBef{j} = [curvPks.dirRelCenterBef{j},...
                curvPksTmp.dirRelCenterBef];
            curvPks.dirRelCenterAft{j} = [curvPks.dirRelCenterAft{j},...
                curvPksTmp.dirRelCenterAft];
            
            % assign to output matrix
            curvWalks.max{j} = [curvWalks.max{j} curvWalksTmp.max];
            curvWalks.tot{j} = [curvWalks.tot{j} curvWalksTmp.tot];
            curvWalks.dur{j} = [curvWalks.dur{j} curvWalksTmp.dur];
            curvWalks.avg{j} = [curvWalks.avg{j} curvWalksTmp.avg];
            curvWalks.ndx{j} = [curvWalks.ndx{j} curvWalksTmp.ndx+startNdx(k)-1];
            curvWalks.all{j} = [curvWalks.all{j} curvWalksTmp.all];
            curvWalks.dirRelCenterBef{j} = [curvWalks.dirRelCenterBef{j},...
                curvWalksTmp.dirRelCenterBef];
            curvWalks.dirRelCenterAft{j} = [curvWalks.dirRelCenterAft{j},...
                curvWalksTmp.dirRelCenterAft];
            
        end
        p = p+2;
    end
    %checkSce(sAll(K,:),curvPks,curvWalks,boundCond,stopCond)
    
    % create an error figure
    figure(figNum);subplot(2,1,1);
    histogram(err);ylabel('Error (% mm/frame)');
    subplot(2,1,2);
    histogram(log10(err));
    xlabel('log Error (% mm/frame)');ylabel('Count')
    suptitle(genAll{i})
    set(gcf,'Position',[2 42 638 600])
    
    % print the figures
    if meta.plotFig
        for f = 1:figNum-1
            figure(f);
            suptitle(genAll{i})
            print('-painters','-dpsc2','allFigures.ps','-loose','-append');
        end
    end
    
    % save the data
    fileName  = [genFolder '\' genAll{i} '_' meta.d '.mat'];
    if isfile(fileName)
         save(fileName,'curvPks','curvWalks','stopCond','boundCond','err','-append','-v7.3');
    else
         save(fileName,'curvPks','curvWalks','stopCond','boundCond','err','-v7.3');
    end
    close all
    errAll{i} = err;
end
%sum(err,2)./sum(err>0,2);
%figure;histogram(log10(ans))

end