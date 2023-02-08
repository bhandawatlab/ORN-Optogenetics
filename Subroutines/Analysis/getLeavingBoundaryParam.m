function [f_orco] = getLeavingBoundaryParam(f_orco,fe,plotFig)
% getLeavingBoundaryParam  calculates the speed, curvature, and angle at
%   which flies leave the physical arena boundary
%
%   Inputs: f_orco = fly object
%           fe = vector of first entry for each fly
%
%   Output: f_orco = fly object updated with the leaving boundary params
%   

%load(['DataModel/' gen '_' meta.d meta.ext '.mat'],'f_orco');

%fe = f_orco.getFirstEntry('H',meta.border);
boundaryNdx = find(cellfun(@(x) strcmpi(x,'boundary'),f_orco.states.key));
spd = f_orco.calcSpd;
avgCurv = cell(f_orco.nFly,1);totCurv = cell(f_orco.nFly,1);
avgSpd = cell(f_orco.nFly,1);dur = cell(f_orco.nFly,1);
leaveStateType = cell(f_orco.nFly,1);rPos = cell(f_orco.nFly,1);
nextState_x = cell(f_orco.nFly,1);nextState_y = cell(f_orco.nFly,1);
nStates2Cons = 3;minTrackLen = 10;
for fly = 1:f_orco.nFly
    currFlyState = f_orco.states.ndx(fly,fe(fly):end);
    currFlyStateTransition = find(1-[1 diff(currFlyState)==0]);
    [~,endNdx,type] = startEndSeq(currFlyState==boundaryNdx);
    %startNdx = startNdx(type==1);
    endNdx = endNdx(type==1);
    endNdx(endNdx==numel(currFlyState)) = [];
    
    
    
    for i = 1:nStates2Cons
        leaveStateType{fly,i} = [currFlyState(endNdx(~isnan(endNdx))+1)'; nan(sum(isnan(endNdx)),1)];
        
        nextStateStart = endNdx+1+fe(fly)-1;% since beginning
        [~,nextStateEnd,idxB] = findBeforeAfter(endNdx+1,currFlyStateTransition,'after');
        nextStateEnd = nextStateEnd-1+fe(fly)-1;
        
        longState = (nextStateEnd-nextStateStart+1)>=minTrackLen;
        
        for j = 1:numel(endNdx)
            
            if longState(j) && ~isnan(nextStateEnd(j)) && nextStateEnd(j)<f_orco.nPt
                avgCurv{fly,i}(j,1) = nanmean(f_orco.curv(fly,nextStateStart(j):nextStateEnd(j))).*180./pi;
                totCurv{fly,i}(j,1) = sum(f_orco.curv(fly,nextStateStart(j):nextStateEnd(j))).*180./pi;
                avgSpd{fly,i}(j,1) = nanmean(spd(fly,nextStateStart(j):nextStateEnd(j)));
                rPos{fly,i}(j,1) = f_orco.rH(fly,nextStateEnd(j));
                nextState_x{fly,j}{i} = f_orco.xH(fly,nextStateStart(j):nextStateEnd(j));
                nextState_y{fly,j}{i} = f_orco.yH(fly,nextStateStart(j):nextStateEnd(j));
            else
                avgCurv{fly,i}(j,1) = nan;
                totCurv{fly,i}(j,1) = nan;
                avgSpd{fly,i}(j,1) = nan;
                rPos{fly,i}(j,1) = nan;
                nextState_x{fly,j}{i} = nan;
                nextState_y{fly,j}{i} = nan;
            end
        end
        
        endNdx = nextStateEnd-fe(fly)+1;
        dur{fly,i} = (nextStateEnd-nextStateStart+1)';
    end
    
    % 1.) leave state type
    % 2.) average speed, curvature, and duration of next state
    % 3.) direction relative to center
    dur{fly,1} = (nextStateEnd-nextStateStart+1)';
end
%%
totCurv_all = cell2mat(totCurv);


tracks2Cons = 1;thresh_step = 15;
figure;
% for t= 1:12
%     subplot(4,3,t)
%     plotCircle([0 0],4,100,'k');hold on;
% end
plotCircle([0 0],4,100,'k');hold on;
leaveAngle = cell(f_orco.nFly,1);
for fly = 1:f_orco.nFly
    for j = 1:find(~cellfun(@isempty,nextState_x(fly,:)), 1, 'last' )
        try
            x_tmp = cell2mat(nextState_x{fly,j}(tracks2Cons));
            y_tmp = cell2mat(nextState_y{fly,j}(tracks2Cons));
            initAngle = myatan(x_tmp(1),y_tmp(1),'degrees',1);
            rotAng = -initAngle;
            
            R = [cosd(rotAng) -sind(rotAng); sind(rotAng) cosd(rotAng)];
            
            xy_rot = R*[x_tmp;y_tmp];
            
            initAng = myatan(diff(xy_rot(1,:)),diff(xy_rot(2,:)),'degrees',1);
            rotAng2 = -nanmean(initAng(1:5));
            R2 = [cosd(rotAng2) -sind(rotAng2); sind(rotAng2) cosd(rotAng2)];
            xy_rot2 = R2*xy_rot;
            xy_rot2 = xy_rot2-xy_rot2(:,1)+xy_rot(:,1);
            
            plot(xy_rot(1,:),xy_rot(2,:));hold on;
            %         for t = 1:12
            %             %tmpCurv = abs(rotAng2);
            %             tmpCurv = abs(totCurv{fly,tracks2Cons}(j));
            %             if tmpCurv>=(thresh_step*(t-1)+1) && tmpCurv<(thresh_step*t)
            %                 subplot(4,3,t)
            %                 plot(xy_rot(1,:),xy_rot(2,:));
            %                 %plot(xy_rot2(1,:),xy_rot2(2,:));
            %             end
            %         end
            leaveAngle{fly,tracks2Cons}(j,1) = rotAng2;
        catch
            a = 1;
            leaveAngle{fly,tracks2Cons}(j,1) = nan;
        end
        
        
    end
end

%%
leaveStateType_all = cell2mat(leaveStateType);
avgCurv_all = cell2mat(avgCurv);
totCurv_all = cell2mat(totCurv);
avgSpd_all = cell2mat(avgSpd);
dur_all = cell2mat(dur);
rPos_all = cell2mat(rPos);
leaveAngle_all = cell2mat(leaveAngle);

corrNdx = (sum(diff(sign(avgCurv_all),[],2)==0))./...
    sum(~(isnan(avgCurv_all(:,1:end-1)) | isnan(avgCurv_all(:,2:end))));

state = 2;
tmp = leaveStateType_all(:,1)==state;

paramsAll = abs([avgCurv_all(tmp,1) avgSpd_all(tmp,1) dur_all(tmp,1)./f_orco.fs leaveAngle_all(tmp,1)]);
paramsAll(any(isnan(paramsAll),2),:) = [];
paramsLab = {'avg curv','avg spd','dur','leave angle'};

figure;k = 1;
for i = 1:size(paramsAll,2)-1
    for j = i+1:size(paramsAll,2)
        subplot(3,2,k);
        scatter(paramsAll(:,i),paramsAll(:,j))
        xlabel(paramsLab{i});ylabel(paramsLab{j})
        title(corr(paramsAll(:,i),paramsAll(:,j), 'rows','complete'))
        k = k+1;
    end
end

%%
% 1.) leave angle
% 2.) average speed, curvature, and duration of next state
figure;
pd = cell(4,1);xyz = zeros(size(paramsAll,2),10);
for param = 1:size(paramsAll,2)
    xx = linspace(min(paramsAll(:,param)),max(paramsAll(:,param)),20);
    pd{param} = fitdist(paramsAll(:,param),'kernel');
    y = pdf(pd{param},xx);
    y_bar = histcounts(paramsAll(:,param),xx,'Normalization','pdf');
    
    subplot(2,2,param);
    bar((xx(1:end-1)+xx(2:end))./2,y_bar);hold on;
    plot(xx,y);
    title(paramsLab{param})
end
f_orco.model.boundary.leaving.avgCurv = pd{1};
f_orco.model.boundary.leaving.spd = pd{2};
f_orco.model.boundary.leaving.dur = pd{3};
f_orco.model.boundary.leaving.angle = pd{4};
f_orco.model.boundary.leaving.corrNdx = corrNdx;

if ~plotFig
    close all
end
%save(['DataModel/' gen '_' meta.d meta.ext '.mat'],'f_orco');
%save('C:\Users\lt532\Desktop\ORN Optogenetics All Code\ORN-Optogenetics-mainFinalJune2022\DataModel\Orco Retinal_April2022_allTime.mat','f_orco')
%save('C:\Users\lt532\Desktop\ORN Optogenetics All Code\ORN-Optogenetics-mainFinalJune2022\DataModel\Max Attraction Retinal_April2022_allTime.mat', 'f_orco')

end
