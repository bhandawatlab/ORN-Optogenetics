classdef Flies < handle
    properties
        id
        nFly
        nPt
        fs
        rBound
        x
        y
        r
        xH
        yH
        rH
        curv
        spk
        lightOn
        states
        model
    end
    methods
        function [self] = Flies(id,x,y,xH,yH,spk,fs,rBound,lightOn,...
                curvPks,curvWalks,stopCond,boundCond,nPt)
            self.id = id;
            self.nFly = size(x,1);
            
            % To account for genotypes with different number of data points
            % tracked, set the length to the minimum length
            missingData = sum(x==0 & y==0 & xH==0 & yH==0)>0;
            if any(missingData) && isempty(nPt)
                self.nPt = find(missingData,1)-1;
            else
                self.nPt = size(x,2);
            end
            self.fs = fs;
            self.rBound = rBound;
            self.x = x(:,1:self.nPt);self.y = y(:,1:self.nPt);
            self.xH = xH(:,1:self.nPt);self.yH = yH(:,1:self.nPt);
            self.r = sqrt(self.x.^2+self.y.^2);
            self.rH = sqrt((self.xH).^2+(self.yH).^2);
            [self.curv,~] = self.calcCurv;
            self.spk = spk(:,1:self.nPt);
            if ~(isempty(curvPks) || isempty(curvWalks) || isempty(stopCond)...
                    || isempty(boundCond))
                [self.states.ndx, self.states.key] = getStates(curvPks,curvWalks,stopCond,boundCond,self.nPt-1);
%                 spd = self.calcSpd;
%                 self.states.ndx(spd<0.5 & self.states.ndx <=2) = 0;
            end
            if isempty(lightOn)
                self.lightOn = ceil(size(x,2)./2).*ones(self.nFly,1);
            else
                self.lightOn = reshape(lightOn,[],1);
            end
            
        end
        
        % calculate curvature
        function [curv,thetaAll] = calcCurv(self)
            curv = zeros(self.nFly,self.nPt-1);
            thetaAll = zeros(self.nFly,self.nPt-1);
            for i = 1:self.nFly
                % calculating curvature
                Vertices = horzcat(self.x(i,:)',self.y(i,:)');
                %get normal vectors
                N=LineNormals2D(Vertices);
                
%                 figure,
%                 plot([Vertices(3720:3750,1) Vertices(3720:3750,1)+N(3720:3750,1)]',[Vertices(3720:3750,2) ...
%                     Vertices(3720:3750,2)+N(3720:3750,2)]');hold on
%                 plot3(Vertices(3720:3750,1),Vertices(3720:3750,2),3720:3750,'k')
%                 view(2)
                
                theta = zeros(1,self.nPt-1); %make it zero vector
                for p=1:(self.nPt-1)
                    %because normal is a unit vector. The x-component can be used
                    % to determine the angle it makes with x-axis.
                    theta(p)=acos(N(p,1));
                    % the if loop converts from 0 to pi to -pi to pi
                    if N(p,1)>0 && N(p,2)<0
                        theta(p)=-theta(p);
                    elseif N(p,1)<0 && N(p,2)<0
                        theta(p)=-theta(p);
                    end
                end
                theta(isnan(theta)) = 0;
                theta = unwrap(theta);
                thetaAll(i,1:end) = theta;
                curv(i,1:end-1) = diff(theta);
            end
            curv = [curv, zeros(self.nFly,1)];
        end
        
        % calculate the speed
        function spd = calcSpd(self)
            spd = sqrt(diff(self.x.*10,[],2).^2+diff(self.y.*10,[],2).^2).*self.fs;
            spd = [spd spd(:,end)];
        end
        
       % calculate the change in firing rate
        function dF = calcDeltaFR(self)
            % df
            dF = gradient(self.spk).*self.fs;
            % df Smooth
        end
        
        % calculate border angle
        function phi = calcPhi(self)
            theta = zeros(size(self.x));
            for fly = 1:self.nFly
                theta(fly,:) = myatan(self.x(fly,:),self.y(fly,:),'degrees',2);
            end
            phi = diff(theta,[],2);
            phi2 = cat(3,abs(phi), abs(abs(phi)-360));
            [~,ndx] = min(phi2,[],3);
            phi(ndx==2) = sign(phi(ndx==2)).*abs(abs(phi(ndx==2))-360);
            phi = [phi phi(:,end)];
        end
        
        % calculate change in orientation
        function [ang,dAng] = calcAng(self)
            ang = zeros(size(self.x));
            for fly = 1:self.nFly
                ang(fly,:) = myatan(self.xH(fly,:)-self.x(fly,:),...
                    self.yH(fly,:)-self.y(fly,:),'degrees',2);
            end
            dAng = diff(ang,[],2);
            dAng2 = cat(3,abs(dAng), abs(abs(dAng)-360));
            [~,ndx] = min(dAng2,[],3);
            dAng(ndx==2) = sign(dAng(ndx==2)).*abs(abs(dAng(ndx==2))-360);
            dAng = [dAng dAng(:,end)];
        end
        
        % get indices for both the condition (what the fly is experiencing)
        % and the movement state
        function [condState,keyNdx] = getCond(self,condNdx,key)
            if any(size(self.states.ndx)~=size(condNdx))
                error('dimensions do not match')
            else
                condState = nan(size(self.states.ndx));
                for i = 1:self.nFly
                    [startNdx,endNdx,type] = startEndSeq(self.states.ndx(i,:));% movement state
%                     [startNdx2,endNdx2,type2] = startEndSeq(condNdx(i,:));% condition spk rate
                    for j = 1:numel(type)
                        [GC,GR] = groupcounts(condNdx(i,startNdx(j):endNdx(j))');
                        [~,tmpCondNdx] = max(GC,[],1);
                        condState(i,startNdx(j):endNdx(j)) = (type(j)-1).*numel(key)+GR(tmpCondNdx);
                    end
                end
                
                %condState = (self.states.ndx-1).*numel(key)+condNdx;
                [A,B] = meshgrid(self.states.key,key);
                c=cat(2,A,B);
                keyNdx=reshape(c,[],2);
            end
        end
        
        function [durState] = getCondDur(self)
            durState = zeros(size(self.states.ndx));
            for fly = 1:self.nFly
                [startNdx,endNdx,type] = startEndSeq(self.states.ndx(fly,:));
                dur = endNdx-startNdx+1;
                for i = 1:numel(startNdx)
                    durState(fly,startNdx(i):endNdx(i)) = dur(i); 
                end
            end
        end
        
        % get the speed, curvature, and time for each condition
        function [s,keyNdx] = getKinematicsCond(self,condNdx,key,thresh,varargin)
            if ~isempty(varargin)
                ndxOnly = varargin{1};
            else
                ndxOnly = [];
            end
            [s,keyNdx] = calcKinematicsCond(self,condNdx,key,thresh,ndxOnly);
        end
        
        % get the first entry
        function fe = getFirstEntry(self,lab,border)
            if isempty(lab)
                rPos = self.rH;
            else
                rPos = self.(['r' lab]);
            end
            in = rPos<border;
            fe = zeros(self.nFly,1);
            for j = 1:self.nFly
                duringIn = in(j,:);
                duringIn(1:self.lightOn(j)-1) = false;
                try
                    fe(j) = find(duringIn,1,'first');
                catch
                    fe(j) = nan;
                end
            end
        end
        
        % remove flies
        function fly2 = rmvData(self,badNdx)
            goodNdx = ~badNdx;
            fly2 = Flies(self.id,self.x(goodNdx,:),self.y(goodNdx,:),...
                self.xH(goodNdx,:),self.yH(goodNdx,:),self.spk(goodNdx,:),...
                self.fs,self.rBound,self.lightOn(goodNdx,:),[],[],[],[],self.nPt);
            if ~isempty(self.states)
                try
                    fly2.states.ndx = self.states.ndx(goodNdx,:);
                catch
                    fly2.states.ndx = self.states.ndx;
                end
                fly2.states.key = self.states.key;
            end
        end
        
        function cross = getCrossing(self,lab,border,minDur,spd)
            % exit/enter
            if numel(border)<2
                border = [border border];
            end
            if isempty(minDur)
                minDur = self.fs;
            end
            t = self.fs;%30;
            fe = self.getFirstEntry('H',border(2));
            if isempty(spd)
                spd = self.calcSpd;
            end
            
            pos = self.(['r' lab]);
            crossType = zeros(self.nFly,self.nPt);
            %crossType(1,pos>border(1)) = 1;%exit
            crossType(pos<border(2)) = 1;%enter
            
            %cross = cell(self.nFly,1);
            for fly = 1:self.nFly
                %enter
                currTrack = crossType(fly,fe(fly):end);
                [startNdx,endNdx,type] = startEndSeq(currTrack);
                % locations after leaving 1 mm
                p = find(pos(fly,fe(fly):end)<border(1));
                
                % find the nearest position after entering the XX cm mark
                % that the fly leaves at the YY cm mark
                try
                [closestPtBef,closestPtAft,idxB] = findBeforeAfter(endNdx,p,'before');
                catch
                    a;
                end
                [currExit,ndx] = unique(closestPtBef);
                % get the absolute time points of the enter and exit
                currExit = reshape(currExit,[],1)+fe(fly)-1;
                currEnter = reshape(startNdx(ndx),[],1)+fe(fly)-1;
                
                currEnter(isnan(currExit)) = [];
                currExit(isnan(currExit)) = [];
                
%                 % plot to verify
%                 figure;plot(fe(fly):self.nPt,pos(fly,fe(fly):end),'k');hold on;
%                 for i = 1:numel(currEnter)
%                 	plot([currEnter(i):currExit(i)],pos(fly,[currEnter(i):currExit(i)]))
%                 end
                
                
                    
                % remove small transitions
                len = currExit-currEnter+1;
                currEnter(len<minDur) = [];currExit(len<minDur) = [];
                
                if ~isempty(currEnter)     
                    cross.ndx{fly} = [currEnter currExit];
                    
                    spkTmp = [self.spk(fly,:) nan(1,t)];
                    spdTmp = [spd(fly,:) nan(1,t)];
                    curvTmp = [self.curv(fly,:) nan(1,t)];
                    posTmp = [pos(fly,:) nan(1,t)];
                    stateTmp = [self.states.ndx(fly,:) nan(1,t)];
                    
                    cross.enter.spk{fly} = spkTmp(currEnter + (-t:t));
                    cross.enter.spd{fly} = spdTmp(currEnter + (-t:t));
                    cross.enter.curv{fly} = curvTmp(currEnter + (-t:t));
                    cross.enter.pos{fly} = posTmp(currEnter + (-t:t));
                    cross.enter.state{fly} = stateTmp(currEnter + (-t:t));
                    
                    cross.exit.spk{fly} = spkTmp(currExit + (-t:t));
                    cross.exit.spd{fly} = spdTmp(currExit + (-t:t));
                    cross.exit.curv{fly} = curvTmp(currExit + (-t:t));
                    cross.exit.pos{fly} = posTmp(currExit + (-t:t));
                    cross.exit.state{fly} = stateTmp(currExit + (-t:t));
                    
                    t2 = t;%30;
                    for i = 1:numel(currEnter)
                        cross.track{fly,i} = [spkTmp(currEnter(i):currExit(i)+t2);...
                            posTmp(currEnter(i):currExit(i)+t2);
                            (1:currExit(i)-currEnter(i)+1+t2);
                            currEnter(i):currExit(i)+t2];
                    end
                    for i = 1:numel(currExit)
                        try
                        cross.trackExit(fly,i) = currEnter(i+1)-currExit(i)+1;
                        catch
                            cross.trackExit(fly,i) = self.nPt-currExit(i)+1;
                        end
                    end
                end
                
            end
            
        end
        
        function cross = getCrossing_before(self,lab,border,minDur,spd)
            % exit/enter
            if numel(border)<2
                border = [border border];
            end
            if isempty(minDur)
                minDur = self.fs;
            end
            t = self.fs;%30;
            fe = self.getFirstEntry('H',border(2));
            if isempty(spd)
                spd = self.calcSpd;
            end
            
            pos = self.(['r' lab]);
            crossType = zeros(self.nFly,self.nPt);
            %crossType(1,pos>border(1)) = 1;%exit
            crossType(pos<border(2)) = 1;%enter
            
            %cross = cell(self.nFly,1);
            for fly = 1:self.nFly
                %enter 
                currTrack = crossType(fly,1:fe(fly));
                [startNdx,endNdx,type] = startEndSeq(currTrack);
                % locations after leaving 1 mm
                p = find(pos(fly,1:fe(fly))<border(1));
                
                % find the nearest position after entering the XX cm mark
                % that the fly leaves at the YY cm mark
                try
                [closestPtBef,closestPtAft,idxB] = findBeforeAfter(endNdx,p,'before');
                catch
                    a;
                end
                [currExit,ndx] = unique(closestPtBef);
                % get the absolute time points of the enter and exit
                currExit = reshape(currExit,[],1);
                currEnter = reshape(startNdx(ndx),[],1);
                
                currEnter(isnan(currExit)) = [];
                currExit(isnan(currExit)) = [];
                
%                 % plot to verify
%                 figure;plot(1:fe(fly),pos(fly,1:fe(fly)),'k');hold on;
%                 for i = 1:numel(currEnter)
%                 	plot([currEnter(i):currExit(i)],pos(fly,[currEnter(i):currExit(i)]))
%                 end
                
                % remove small transitions
                len = currExit-currEnter+1;
                currEnter(len<minDur) = [];currExit(len<minDur) = [];
                
                if ~isempty(currEnter)
                    cross.ndx{fly} = [currEnter currExit];

                    spkTmp = [self.spk(fly,:) nan(1,t)];
                    spdTmp = [spd(fly,:) nan(1,t)];
                    curvTmp = [self.curv(fly,:) nan(1,t)];
                    posTmp = [pos(fly,:) nan(1,t)];
                    stateTmp = [self.states.ndx(fly,:) nan(1,t)];

                    enterNdx = currEnter + (-t:t);
                    exitNdx = currExit + (-t:t);
                    cross.enter.spk{fly} = spkTmp(max(enterNdx,1));
                    cross.enter.spd{fly} = spdTmp(max(enterNdx,1));
                    cross.enter.curv{fly} = curvTmp(max(enterNdx,1));
                    cross.enter.pos{fly} = posTmp(max(enterNdx,1));
                    cross.enter.state{fly} = stateTmp(max(enterNdx,1));
                    
                    cross.exit.spk{fly} = spkTmp(max(exitNdx,1));
                    cross.exit.spd{fly} = spdTmp(max(exitNdx,1));
                    cross.exit.curv{fly} = curvTmp(max(exitNdx,1));
                    cross.exit.pos{fly} = posTmp(max(exitNdx,1));
                    cross.exit.state{fly} = stateTmp(max(exitNdx,1));
                    
                    cross.enter.spk{fly}(enterNdx<1) = nan;
                    cross.enter.spd{fly}(enterNdx<1) = nan;
                    cross.enter.curv{fly}(enterNdx<1) = nan;
                    cross.enter.pos{fly}(enterNdx<1) = nan;
                    cross.enter.state{fly}(enterNdx<1) = nan;
                    
                    cross.exit.spk{fly}(exitNdx<1) = nan;
                    cross.exit.spd{fly}(exitNdx<1) = nan;
                    cross.exit.curv{fly}(exitNdx<1) = nan;
                    cross.exit.pos{fly}(exitNdx<1) = nan;
                    cross.exit.state{fly}(exitNdx<1) = nan;
                    
                    t2 = t;%30;
                    for i = 1:numel(currEnter)
                        cross.track{fly,i} = [spkTmp(currEnter(i):currExit(i)+t2);...
                            posTmp(currEnter(i):currExit(i)+t2);
                            (1:currExit(i)-currEnter(i)+1+t2);
                            currEnter(i):currExit(i)+t2];
                    end
                end
            end
            
        end
        
        function cross = getCrossingBaseline(self,lab,border,baseline,minDur)
            fe = self.getFirstEntry('H',border);
            spd = self.calcSpd;
            pos = self.(['r' lab]);
            baseLineNdx = (self.spk-baseline)<0.001;%abs(self.spk-baseline)<0.001;
            t = self.fs;%30
            
            if isempty(minDur)
                minDur = self.fs;
            end
            
            for fly = 1:self.nFly
                currTrack = baseLineNdx(fly,fe(fly):end);
                [startNdx,endNdx,type] = startEndSeq(currTrack);
                startNdx = startNdx(type==0);
                endNdx = endNdx(type==0);
                len = endNdx-startNdx+1;
                startNdx(len<minDur) = [];endNdx(len<minDur) = [];
                
                currExit = reshape(endNdx,[],1)+fe(fly)-1;
                currEnter = reshape(startNdx,[],1)+fe(fly)-1;
                
                spkTmp = [self.spk(fly,:) nan(1,t)];
                spdTmp = [spd(fly,:) nan(1,t)];
                curvTmp = [self.curv(fly,:) nan(1,t)];
                posTmp = [pos(fly,:) nan(1,t)];
                stateTmp = [self.states.ndx(fly,:) nan(1,t)];
                
                cross.enter.spk{fly} = spkTmp(currEnter + (-t:t));
                cross.enter.spd{fly} = spdTmp(currEnter + (-t:t));
                cross.enter.curv{fly} = curvTmp(currEnter + (-t:t));
                cross.enter.pos{fly} = posTmp(currEnter + (-t:t));
                cross.enter.state{fly} = stateTmp(currEnter + (-t:t));
                
                cross.exit.spk{fly} = spkTmp(currExit + (-t:t));
                cross.exit.spd{fly} = spdTmp(currExit + (-t:t));
                cross.exit.pos{fly} = posTmp(currExit + (-t:t));
                cross.exit.curv{fly} = curvTmp(currExit + (-t:t));
                cross.exit.state{fly} = stateTmp(currExit + (-t:t));
                
                for i = 1:numel(currEnter)
%                     if i == 3 && fly == 4
%                         figure;plot(spkTmp(currEnter(i)-30:currExit(i)))
%                     end
                    cross.track{fly,i} = [spkTmp(currEnter(i):currExit(i));...
                        posTmp(currEnter(i):currExit(i));
                        (1:currExit(i)-currEnter(i)+1);
                        currEnter(i):currExit(i)];
                end
            end
            
        end
        
        % get attraction index
        function attractNdx = getAttnNdx(self,lab,border)
            attractNdx = calcAttractionIndex(self,lab,border);
        end
        
        % get radialProbability
        function radialProb = getRadialProb(self,lab,border,stopSpd,dur)
            radialProb = calcRadialProb(self,lab,border,stopSpd,dur);
        end
        
        % get the turn density
        function turnDens = getTurnDens(self,lab,s,key)
            turnDens = calcTurnDens(self,lab,s,key);
        end
        
        % get the prob of being inside
        function [ProbIn] = getProbIn(self,lab,border,startPts)
            [bin_inside_rim_sliding_SEM,bin_inside_rim_sliding_mean,inside_rim,framesinBin]...
                = calcProbIn(self,lab,border,startPts);
            ProbIn.sem = bin_inside_rim_sliding_SEM;
            ProbIn.mean = bin_inside_rim_sliding_mean;
            ProbIn.inside_rim = inside_rim;
            ProbIn.framesinBin = framesinBin;
        end
        
        % get state kinematics and turn optimality based on ORN activity
        function [] = GetKinematicModelParams(self,condNdx,key,fe,opts,plotFig)
            % KNN distribution for when there are dynamic ORN activity
            getDistByTimeSpike(self,condNdx,key,opts,plotFig);
            % Distribution when ORN activity reaches baseline again
            getBaselineKinematicsEarlyLate(self,opts.border,condNdx,key);%getBaselineKinematics(self,condNdx,key)
            % arena boundary condition
            getMeanBound(self,fe);
        end
        
    end
end