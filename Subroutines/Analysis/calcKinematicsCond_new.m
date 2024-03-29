function [s,keyNdx] = calcKinematicsCond_new(self,condNdx,key,thresh,ndxOnly,timeInterval)
% calcKinematicsCond_new  divides each state instances into cell arrays
%   based on firing rate regime. Key descriptors of each state such as
%   kinematics, index, prior firing rate info, and turn optimality are
%   outputed in a structure.
%
%   Inputs: self = fly object
%           condNdx = fly x time matric of the general firing rate regime
%           the fly is in (1 = non-baseline, 2 = before, 3 = during baseline)
%           key = cell array key for condNdx
%           thresh = percent threshold classified into main condition to
%               consider (i.e. if 0, then only tracks that are 100 % in a
%               condition is considered for that condition)
%           ndxOnly = true/false to calculate all information about states
%               in each condition or only return the index of each state
%           timeInterval = time period to average over before start of each
%               state trajectory (in ms)
%
%   Output: s = structure of descriptions for each state (e.g. average 
%               curvature, average speed, time ndx, prior ORN firing rate, 
%               etc). Each field (description) is a fly x condition cell
%               array. Each cell is a m x 1 vector of descriptors for each
%               state instance
%           keyNdx = Cell array describing the condition (each column) of
%               the cell arrays
%   

nKey = numel(key);nStates = numel(self.states.key);
bins = 1:nKey+1;
spd = self.calcSpd;
[~, thetaAll] = self.calcCurv;
thetaAll = thetaAll+pi;
[ang,dAng] = self.calcAng;

[s] = initS(self.nFly,nStates*nKey);
dSpkRateAll = self.calcDeltaFR;
dt = diff(timeInterval);
[dFSmooth_dt,FSmooth_dt,~] = get_df_f_history(self,dSpkRateAll,[0 dt]);
offSet = round(timeInterval(1)./1000.*self.fs);
dFSmooth = [zeros(self.nFly,offSet),dFSmooth_dt(:,1:end-offSet)];
FSmooth = [zeros(self.nFly,offSet),FSmooth_dt(:,1:end-offSet)];

% [dFSmooth_test,FSmooth_test,~] = get_df_f2(self,dSpkRateAll,timeInterval);
% assert(all((dFSmooth_test-dFSmooth)==0,'all'))
% assert(all((FSmooth_test-FSmooth)==0,'all'))

hist = ceil(0.2.*self.fs)-1;%200 ms history

for fly = 1:self.nFly
    [startNdx,endNdx,type] = startEndSeq(self.states.ndx(fly,:));
    startNdx(type == 0) = [];endNdx(type == 0) = [];type(type == 0) = [];
    startNdx = reshape(startNdx,[],1);endNdx = reshape(endNdx,[],1);
    for i = 1:nStates
        startNdxTmp = startNdx(type == i);
        endNdxTmp = endNdx(type == i);
        len = endNdxTmp-startNdxTmp+1;
        percContr = zeros(numel(startNdxTmp),nKey);
        for j = 1:numel(startNdxTmp)
            percContr(j,:) = histcounts(condNdx(fly,startNdxTmp(j):endNdxTmp(j)),bins)./len(j);
        end
        badTracks = ~any(percContr>thresh,2) | len<3;
        startNdxTmp(badTracks) = [];
        endNdxTmp(badTracks) = [];
        percContr(badTracks,:) = [];
        [~,c] = max(percContr,[],2);
        
        %%%
%         badTracks_lastBaseline = false(numel(startNdxTmp),1);
%         for j = 1:numel(startNdxTmp)
%             if c(j) ==1 % not before first entry
%                 F_hist = FSmooth_dt(fly,startNdxTmp(j)-offSet:startNdxTmp(j));
%                 dF_hist = dFSmooth_dt(fly,startNdxTmp(j)-offSet:startNdxTmp(j));
%                 badTracks_lastBaseline(j) = ...
%                     (any(abs(F_hist-self.spk(1))<0.01) || any(abs(dF_hist)<0.01));
%             end
%         end
%         
%         startNdxTmp(badTracks_lastBaseline) = [];
%         endNdxTmp(badTracks_lastBaseline) = [];
%         percContr(badTracks_lastBaseline,:) = [];
%         [~,c] = max(percContr,[],2);
        %%%
        
        kk = ones(1,nStates.*nKey);
        if ndxOnly == true
            for j = 1:numel(startNdxTmp)
                k = (i-1).*nKey+c(j);
                s.ndx{fly,k}(kk(k),1) = startNdxTmp(j)+1;
                kk(k) = kk(k)+1;
            end
        else
            for j = 1:numel(startNdxTmp)
                k = (i-1).*nKey+c(j);

                spdTmp = spd(fly,startNdxTmp(j):endNdxTmp(j));
                s.avgSpd{fly,k}(kk(k),1) = nanmean(spdTmp);
                s.totSpd{fly,k}(kk(k),1) = nansum(spdTmp);
                s.stdSpd{fly,k}(kk(k),1) = nanstd(spdTmp);
                curvTmp = self.curv(fly,startNdxTmp(j):endNdxTmp(j));
                s.avgCurv{fly,k}(kk(k),1) = nanmean(curvTmp);
                s.totCurv{fly,k}(kk(k),1) = nansum(curvTmp);
                s.maxCurv{fly,k}(kk(k),1) = max(abs(curvTmp));

                s.diffSpd{fly,k}(kk(k),1) = (spdTmp(end)-spdTmp(1))./numel(spdTmp);
                s.diffCurv{fly,k}(kk(k),1) = (curvTmp(end)-curvTmp(1))./numel(spdTmp);

                % previous state
                if startNdxTmp(j) == startNdx(1)
                    prevStartNdx = max(startNdx(1),1);
                    prevEndNdx = prevStartNdx+3;

                    s.prevState{fly,k}(kk(k),1) = nan;
                    s.currState{fly,k}(kk(k),1) = nan;
                else
                    prevIndex = find(startNdx==startNdxTmp(j))-1;% previous state
                    prevStartNdx = startNdx(prevIndex);
                    prevEndNdx = startNdxTmp(j)-1;

                    s.prevState{fly,k}(kk(k),1) = type(prevIndex);
                    s.currState{fly,k}(kk(k),1) = type(prevIndex+1);
                end

                % next state
                if endNdxTmp(j) == endNdx(end)
                    s.nextState{fly,k}(kk(k),1) = nan;
                else
                    nextIndex = find(endNdx==endNdxTmp(j))+1;% previous state
                    s.nextState{fly,k}(kk(k),1) = type(nextIndex);
                end

                curvTmp_prev = self.curv(fly,prevStartNdx:prevEndNdx);
                s.avgCurv_prev{fly,k}(kk(k),1) = mean(curvTmp_prev);
                s.totCurv_prev{fly,k}(kk(k),1) = sum(curvTmp_prev);

                %%%%
                startAng = ang(fly,startNdxTmp(j));
                endAng = ang(fly,endNdxTmp(j));
                startDir = [cosd(startAng);sind(startAng)];
                endDir = [cosd(endAng);sind(endAng)];
                startDirRelCent = -[self.x(fly,startNdxTmp(j)); self.y(fly,startNdxTmp(j))];%xH,yH
                endDirRelCent = -[self.x(fly,endNdxTmp(j)); self.y(fly,endNdxTmp(j))];

                n = [0;0;1];
                v = [startDir;0];
                u = [startDirRelCent;0];
                thetaStart = get2VecAngle(v,u,n,false);
                n = [0;0;1];
                v = [endDir;0];
                u = [endDirRelCent;0];
                thetaEnd = get2VecAngle(v,u,n,false);
                thetaDiff = thetaEnd-thetaStart;% negative means turned more inwards
                s.optimalStop{fly,k}(kk(k),1) = thetaDiff;
                %%%%

                % calculate current direction of movement relative to odor
                % source (+ is to the right, - is to the left)
                % if -, right is optimal dir; if +, left is optimal dir
                currAng = thetaAll(fly,max(prevEndNdx-hist,1):prevEndNdx);
                currAng = nanmean(currAng);
                currAng = wrapTo360(currAng.*180./pi)-90;
                prevDir = [cosd(currAng);sind(currAng)];
                dirRelCent = -[self.x(fly,prevEndNdx); self.y(fly,prevEndNdx)];
                n = [0;0;1];
                v = [prevDir;0];
                u = [dirRelCent;0];
                theta = get2VecAngle(v,u,n,false);

                % calculate direction of movement relative to odor source after
                % turning (+ is to the right, - is to the left)
                % if -, right is optimal dir; if +, left is optimal dir
                nextAng = thetaAll(fly,endNdxTmp(j):min(endNdxTmp(j)+1,self.nPt-1));
                nextAng = nanmean(nextAng);
                nextAng = wrapTo360(nextAng.*180./pi)-90;
                nextDir = [cosd(nextAng);sind(nextAng)];
                dirRelCent = -[self.x(fly,prevEndNdx); self.y(fly,prevEndNdx)];
                n = [0;0;1];
                v = [nextDir;0];
                u = [dirRelCent;0];
                theta2 = get2VecAngle(v,u,n,false);

                %--------------------------------------------------------------
                % code for turn optimality based on orientation
                %--------------------------------------------------------------
    %             prevDir = [self.xH(fly,max(prevEndNdx-2,1):prevEndNdx)-...
    %                 self.x(fly,max(prevEndNdx-2,1):prevEndNdx) ;...
    %                     self.yH(fly,max(prevEndNdx-2,1):prevEndNdx)-...
    %                     self.y(fly,max(prevEndNdx-2,1):prevEndNdx)];

                prevDir = [self.xH(fly,startNdxTmp(j):startNdxTmp(j)+2)-...
                    self.x(fly,startNdxTmp(j):startNdxTmp(j)+2) ;...
                        self.yH(fly,startNdxTmp(j):startNdxTmp(j)+2)-...
                        self.y(fly,startNdxTmp(j):startNdxTmp(j)+2)];    

                prevDir = nanmean(prevDir,2);
                dirRelCent = -[self.x(fly,startNdxTmp(j)); self.y(fly,startNdxTmp(j))];
                %dirRelCent = -[self.x(fly,prevEndNdx); self.y(fly,prevEndNdx)];
                n = [0;0;1];
                v = [prevDir;0];
                u = [dirRelCent;0];
                theta_orient = get2VecAngle(v,u,n,false);
                %--------------------------------------------------------------
                theta3 = abs(theta);
                if theta3>180
                    pause(0.1);
                end
                nonOpt = sign(nanmean(theta)) ~= sign(s.avgCurv{fly,k}(kk(k),1));
                theta3(nonOpt) = 360-theta3(nonOpt);
                s.theta{fly,k}(kk(k),1) = theta3;


                % curvature (-) = right turn, (+) = left turn ==> optimal is
                % where curvature sign matches the current direction of
                % movement relative to odor
                s.optimal{fly,k}(kk(k),1) = sign(nanmean(theta)) == sign(s.avgCurv{fly,k}(kk(k),1));
                s.turnBias{fly,k}(kk(k),1) = abs(theta2)<abs(theta);

                s.optimalPLOS{fly,k}(kk(k),1) = sign(nanmean(theta_orient)) == sign(s.avgCurv{fly,k}(kk(k),1));

                s.spkStart{fly,k}(kk(k),1) = FSmooth(fly,startNdxTmp(j));
                %s.spkMid{fly,k}(kk(k),1) = mean(spkRate(midIdx-1:min(midIdx+1,length(spkRate))));
                s.dSpkStart{fly,k}(kk(k),1) = dFSmooth(fly,startNdxTmp(j));
                %s.dSpkMid{fly,k}(kk(k),1) = mean(dSpkRate(midIdx-1:min(midIdx+1,length(spkRate))));
                s.spkEnd{fly,k}(kk(k),1) = FSmooth(fly,endNdxTmp(j));
                s.dSpkEnd{fly,k}(kk(k),1) = dFSmooth(fly,endNdxTmp(j));

    %             if any(spkRate(1:end-hist)==0)
    %                 prevTrack = max(prevEndNdx-hist,1):prevEndNdx;
    %                 figure(100);
    %                 plot(self.x(fly,prevTrack),self.y(fly,prevTrack),'b','Linewidth',2);hold on;
    %                 plot(self.x(fly,startNdxTmp(j):endNdxTmp(j)),self.y(fly,startNdxTmp(j):endNdxTmp(j)),'r','Linewidth',2);
    %                 q0 = quiver(self.x(fly,prevTrack(end)),self.y(fly,prevTrack(end)),-self.x(fly,prevTrack(end)),-self.y(fly,prevTrack(end)));
    %                 q0.Color = 'blue';q0.MaxHeadSize = 1;q0.LineStyle = '--';
    %                 q1 = quiver(self.x(fly,prevTrack(end)),self.y(fly,prevTrack(end)),cosd(currAng)./2,sind(currAng)./2);
    %                 q2 = quiver(self.x(fly,endNdxTmp(j)),self.y(fly,endNdxTmp(j)),cosd(nextAng)./2,sind(nextAng)./2);
    %                 q1.Color = 'blue';q1.MaxHeadSize = 1;q1.LineStyle = '--';
    %                 q2.Color = 'red';q2.MaxHeadSize = 1;q2.LineStyle = '--';
    %                 plotCircle([0,0],1,100,'r');
    %                 plotCircle([0,0],0.01,100,'r');
    %                 axis([-2 2 -2 2]);axis square
    %                 hold off;
    %                 title({['turnBias: ' num2str(s.turnBias{fly,k}(kk(k),1))...
    %                     ' optimality: ' num2str(s.optimal{fly,k}(kk(k),1))...
    %                     ' curv: ' num2str(s.totCurv{fly,k}(kk(k),1))...
    %                     ' theta: ' num2str(s.theta{fly,k}(kk(k),1))], ...
    %                     [' spkBefore: ' num2str(s.spkStart{fly,k}(kk(k),1))...
    %                     ' spkAft: ' num2str(s.spkEnd{fly,k}(kk(k),1))]})
    %                 
    %                 pause(0.1);
    %             end



                try
                    s.locStart{fly,k}(kk(k),1) = mean(self.rH(fly,startNdxTmp(j):startNdxTmp(j)+1));%+l
                catch
                    a;
                end
                                            %mod(abs(sum(abs(curvTmp))).*180./pi,360);
                %------------------
                %s.ndx{fly,k}(kk(k),1) = startNdxTmp(j)+find(abs(curvTmp)==max(abs(curvTmp)),1)-1;
                s.ndx{fly,k}(kk(k),1) = startNdxTmp(j)+1;
                %------------------
                s.dur{fly,k}(kk(k),1) = endNdxTmp(j)-startNdxTmp(j)+1;
                kk(k) = kk(k)+1;
            end
        end
    end
end
[A,B] = meshgrid(self.states.key,key);
c=cat(2,A,B);
keyNdx=reshape(c,[],2);

% figure;hold on
% for i = 1:2
%     currRLoc = cell2mat(s.locStart(:,i));
%     xx =[0:0.1:4.1];
%     plot(histcounts(currRLoc,xx)./(pi.*(xx(2:end).^2-xx(1:end-1).^2)));
% end
% 
% figure;hold on
% for i = 1:1
%     currRLoc = cell2mat(s.locStart(:,i));
%     xx =[0:0.1:4.1];
%     tmp = histcounts(currRLoc,xx)./(pi.*(xx(2:end).^2-xx(1:end-1).^2));
%     plot(xx(1:end-1),tmp./sum(tmp));
%     xlabel('radial position');ylabel('Turn Density')
% end
% 
% s2 = s;
% xx2 = xx;
% load('tempRadLocTurns_1230', 's');
% 
% figure;hold on
% for i = 1:1
%     currRLoc = cell2mat(s.locStart(:,i));
%     xx =[0:0.1:4.1];
%     tmp = histcounts(currRLoc,xx)./(pi.*(xx(2:end).^2-xx(1:end-1).^2));
%     plot(xx(1:end-1),tmp./sum(tmp));hold on;
%     
%     currRLoc = cell2mat(s2.locStart(:,i));
%     tmp = histcounts(currRLoc,xx2)./(pi.*(xx2(2:end).^2-xx2(1:end-1).^2));
%     plot(xx2(1:end-1),tmp./sum(tmp));
%     
%     xlabel('radial position');ylabel('Turn Density')
% end

end


function [s] = initS(nFly,nCond)

s.avgSpd = cell(nFly,nCond);
s.totSpd = cell(nFly,nCond);
s.avgCurv = cell(nFly,nCond);
s.totCurv = cell(nFly,nCond);
s.ndx = cell(nFly,nCond);
s.dur = cell(nFly,nCond);

end

% function [dFSmooth,FSmooth,ndx] = get_df_f2(self,dSpkRateAll,timeInterval)
% 
% hist_start = ceil(timeInterval(1)./1000.*self.fs);%max(1,
% hist_end = ceil(timeInterval(2)./1000.*self.fs)-1;
% dt = hist_end-hist_start+1;
% 
% ndx = [1:size(dSpkRateAll,2)]-[hist_start:hist_end]'-1;
% ndx(:,sum(ndx<1)>0) = [];
% 
% %ndx = (hist_start:size(dSpkRateAll,2)-dt)+[0:dt]'-hist_end;
% for fly = 1:self.nFly
%     curr_df = dSpkRateAll(fly,:);
%     curr_f = self.spk(fly,:);
%     
%     dFSmooth(fly,:) = mean(curr_df(ndx),1);
%     FSmooth(fly,:) = mean(curr_f(ndx),1);
% end
% dFSmooth = [zeros(size(dFSmooth,1),hist_end+1),dFSmooth];
% FSmooth = [zeros(size(FSmooth,1),hist_end+1),FSmooth];
% 
% % %dFSmooth = [zeros(size(dFSmooth,1),hist_end),dFSmooth]./(hist_end-hist_start+1);% average
% % %FSmooth = [zeros(size(FSmooth,1),hist_end),FSmooth]./(hist_end-hist_start+1);% average
% 
% end

% function [dFSmooth,FSmooth,tmp] = get_df_f(self,dSpkRateAll,timeInterval)
% 
% hist_start = max(1,ceil(timeInterval(1)./1000.*self.fs));
% hist_end = ceil(timeInterval(2)./1000.*self.fs)-1;
% dFSmooth = dSpkRateAll(:,1:end-hist_end);
% FSmooth = self.spk(:,1:end-hist_end);
% k = 1;%tmp = [];
% for i = hist_start:hist_end
%     tmp(k,:) = i:size(dSpkRateAll,2)-(hist_end+1)+i;
%     dFSmooth = dSpkRateAll(:,tmp(k,:))+dFSmooth;
%     FSmooth = self.spk(:,tmp(k,:))+FSmooth;
%     %dFSmooth = dSpkRateAll(:,i:end-(hist_end+1)+i)+dFSmooth;
%     %FSmooth = self.spk(:,i:end-(hist_end+1)+i)+FSmooth;
%     k= k+1;
% end
% %dFSmooth = [zeros(size(dFSmooth,1),hist_end),dFSmooth]./(hist_end-hist_start+1);% average
% %FSmooth = [zeros(size(FSmooth,1),hist_end),FSmooth]./(hist_end-hist_start+1);% average
% 
% end