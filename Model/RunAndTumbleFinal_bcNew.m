function [synthFlys,state,nTurn,spdAll,curvAll,phiAll,dSpk,allAng,lambda_CW,...
    lambda_Stop,dfSmooth,fSmooth] = RunAndTumbleFinal_bcNew(self,params)
% get basic parameters
originalFs = self.fs;
fs = params.LLFfs;
len = params.len*params.LLFfs/params.fs;
totFly = params.totFly;
posInit = params.posInit;
border = params.border;
lightOnTime = params.lightOnTime*params.LLFfs/params.fs;
bLFP = params.bLFP;
bRate = params.bRate;
nTotalFilt = (numel(bLFP)+numel(bRate)-2);
ratio = fs./originalFs;
hist = ceil(0.2.*fs)-1;% 200 ms sampling rate

% border choice
C = params.C;
tau = params.tau;

% initialize transition probabilities
TP = self.model.TP;
if size(TP.during,1)<4
    TP.during = [TP.during;TP.during(3,:)];
end
if size(TP.before,1)<4
    TP.before = [TP.before;TP.before(3,:)];
end
% initialize index of states
allStates = self.states.key;
ST_Ndx = find(cellfun(@(c)strcmp(c,{'sharp turns'}),allStates,'UniformOutput',true));
CW_Ndx = find(cellfun(@(c)strcmp(c,{'curved walks'}),allStates,'UniformOutput',true));
Stop_Ndx = find(cellfun(@(c)strcmp(c,{'stops'}),allStates,'UniformOutput',true));
Bound_Ndx = find(cellfun(@(c)strcmp(c,{'boundary'}),allStates,'UniformOutput',true));

% extract the state and kinematic type from each model parameter cell
stateAll = cell(numel(self.model.params),1);
for i = 1:numel(self.model.params)
    tmpState = self.model.params{i}.state.state;
    tmpKin = self.model.params{i}.state.kin;
    stateAll{i} = [tmpState '_' tmpKin];
end
% get the index of the cell that correspond to each kinematic type
ST_spd = find(cellfun(@(c)strcmp(c,{'sharp turns_spd'}),stateAll,'UniformOutput',true));
CW_spd = find(cellfun(@(c)strcmp(c,{'curved walks_spd'}),stateAll,'UniformOutput',true));
ST_curv = find(cellfun(@(c)strcmp(c,{'sharp turns_totCurv'}),stateAll,'UniformOutput',true));
CW_curv = find(cellfun(@(c)strcmp(c,{'curved walks_avgCurv'}),stateAll,'UniformOutput',true));
ST_dur = find(cellfun(@(c)strcmp(c,{'sharp turns_dur'}),stateAll,'UniformOutput',true));
CW_dur = find(cellfun(@(c)strcmp(c,{'curved walks_dur'}),stateAll,'UniformOutput',true));
Stop_curv = find(cellfun(@(c)strcmp(c,{'stops_totCurv'}),stateAll,'UniformOutput',true));
Stop_dur = find(cellfun(@(c)strcmp(c,{'stops_dur'}),stateAll,'UniformOutput',true));

turnBias_ST = find(cellfun(@(c)strcmp(c,{'sharp turns'}),self.model.TurnBias.key,'UniformOutput',true));
turnBias_CW = find(cellfun(@(c)strcmp(c,{'curved walks'}),self.model.TurnBias.key,'UniformOutput',true));
turnBias_Stop = find(cellfun(@(c)strcmp(c,{'stops'}),self.model.TurnBias.key,'UniformOutput',true));


%assume that the fly always start off either in or near the center
xPos = zeros(totFly,len);
yPos = zeros(totFly,len);
rPos = zeros(totFly,len);
if ~isempty(posInit)
    xPos(:,1) = params.posInit(:,1);
    yPos(:,1) = params.posInit(:,2);
end
rPos(:,1) = sqrt(xPos(:,1).^2+yPos(:,1).^2);

% initialize matrices to keep track of angles, locations, etc. Many of
% these are used for debugging purposes
currAng = zeros(totFly,1); allAng = zeros(totFly,len);
during = false(totFly,1);
currDur = zeros(totFly,1);
oldDur = zeros(totFly,1);
currCurv = zeros(totFly,1);
loc = zeros(totFly,len);
state = zeros(totFly,len);
currSpd = zeros(totFly,1);
currPhi = zeros(totFly,1);
newTrack = true(totFly,1);
duringAll = false(totFly,len);
turn = false(totFly,1);
newTrackState = CW_Ndx*ones(totFly,1);% start in the curve walk state
state(:,1) = newTrackState;
nTurn = zeros(totFly,len);
t_relFE = -ones(totFly,len)./fs;
I = zeros(totFly,len);
currSpk = zeros(totFly,len);
dSpk = zeros(totFly,len);
tt_baseline = -1.*ones(totFly,len);
tt_inhibition = -1.*ones(totFly,len);

ttInh = nan(totFly,1);

tSinceOn = inf(totFly,len);
tSinceOff = inf(totFly,len);


dfSmooth = zeros(totFly,len);
fSmooth = zeros(totFly,len);
lambda_CW = zeros(totFly,len);
lambda_Stop = zeros(totFly,len);

nTurn = zeros(totFly,len);

OffState = zeros(totFly,len);
OnState = zeros(totFly,len);

spdAll = zeros(totFly,len);
curvAll = zeros(totFly,len);
phiAll = zeros(totFly,len);
nanVec = nan(totFly,1);

fBaseline = params.baseLineFR;%4.7014;

% generate synthetic tracks based on run and tumble model
progressbar
for i = 1:len-1
    
    % calculate the experienced intensity and resulting spike rate
    %--------------
    if i>lightOnTime
        I(:,i) = location2Intensity(rPos(:,i),self.model.IntensitySpace.x./4,self.model.IntensitySpace.I,self.model.IntensitySpace.convIV);
    end
    tmpSpk = getSpikeRate(I(:,max(1,i-nTotalFilt):i),totFly,bLFP,bRate);
    currSpk(:,i) = tmpSpk(:,end);
    %--------------
    
    % compute the change in spike rate using gradient
    %--------------
    tmpGradient = gradient(currSpk(:,max(i-2,1):i)).*fs;%ratio
    if i>3
        dSpk(:,i-1:i) = tmpGradient(:,2:3);
    elseif i==2
        dSpk(:,i-1:i) = tmpGradient(:,1:2);
    end
    %--------------
    
    inside = currSpk(:,i)>10;% if spike rate is higher than 10 hz, then the fly is inside
    bound = rPos(:,i)>0.985;
    outside = ~bound & ~inside;
    loc(inside,i) = 2;loc(outside,i) = 3;loc(bound,i) = 1;
    
    %----------------------------------------------------------------------
    % this section makes the fly choose a new track if it reaches the edge
    if i>1
        newLoc = (loc(:,i-1) ~= loc(:,i));
        newTrack(newLoc & loc(:,i)==1) = true;
    end
    
    % set flies to the during case after first entry
    if ~all(during) && i>lightOnTime
        % light turns on halfway through
        during(during==false) = inside(during==false);
        duringAll(:,i) = during;
    end
    before = ~during;
    
    if i>1
        t_relFE(during,i) = t_relFE(during,i-1)+1./fs;
        tt_baseline(:,i) = tt_baseline(:,i-1);
    end
    
    % this section makes the fly choose a new track if the fly is in the
    % light ring when the light suddenly turns on
    if i>lightOnTime-1%1
        newLoc = (duringAll(:,i-1) ~= duringAll(:,i));
        newTrack(newLoc & loc(:,i)==2) = true;
    end
    
    % if the current state ends, then a new track begins
    newTrack(currDur<=0) = true;
    %----------------------------------------------------------------------
    
    % this section is for the turn counter
    if i>1
        nTurn(:,i) = nTurn(:,i-1);              % initialize current turn counter as the previous turn number
        tmp = (loc(:,i)==2 & loc(:,i-1)==3) | (loc(:,i)==3 & loc(:,i-1)==2);
        nTurn(tmp & during,i) = 0;              % reset the turn counter if the fly crosses the light border in the during case
    end
    
    %----------------------------------------------------------------------
    % this section is for choosing what new state to enter
    % turn: 1, run: 2, stop: 3, boundary: 4
    if i>1
        f = nanmean(currSpk(:,max((i-hist-1),1):i-1),2);
        df = nanmean(dSpk(:,max((i-hist-1),1):i-1),2);
        df(isnan(df)) = 0;f(isnan(f)) = currSpk(1);% set nans to default 0 hz change and baseline spiking
        
        dfSmooth(:,i) = df;
        fSmooth(:,i) = f;
        df(abs(df)>150) = sign(df(abs(df)>150)).*150;
        
        CW_States = state(:,i-1)==CW_Ndx;
        
        f = nanmean(currSpk(:,max((i-hist-1),1):i-1),2);
        df = nanmean(dSpk(:,max((i-hist-1),1):i-1),2);
        df(isnan(df)) = 0;f(isnan(f)) = currSpk(1);% set nans to default 0 hz change and baseline spiking
        
        dfSmooth(:,i) = df;
        fSmooth(:,i) = f;
        df(abs(df)>150) = sign(df(abs(df)>150)).*150;

        OffState(:,i) = dfSmooth(:,i)<-15;%BC.Thresh(2);        
        tSinceOff(:,i) = tSinceOff(:,i-1)+1;
        tSinceOff(OffState(:,i)==1,i) = 0;
        
        OnState(:,i) = dfSmooth(:,i)>15;%BC.Thresh(1);        
        tSinceOn(:,i) = tSinceOn(:,i-1)+1;
        tSinceOn(OffState(:,i)==1,i) = 0;
        
        % This section is for curved walk border choice
        %----------------------------------------------------
        % leaving
%         OffRate = C.*exp(-(tSinceOff(:,i)./fs)./tau);
%         r = rand(totFly,1)<(OffRate);%./ratio%.*30./100
%         newTrack(r & CW_States) = true;
        % entering
        OnRate = C.*exp(-(tSinceOn(:,i)./fs)./tau);
        r = rand(totFly,1)<(OnRate);%./ratio%.*30./100
        newTrack(r & CW_States) = true;
        %----------------------------------------------------
        
%         % this section makes the fly have an increase chance of choosing a new 
%         % track right when entering/exiting the odor zone
%         if sum(during)>0
%             pSamp = rand(totFly,1);
%             for j = 1:numel(BC_x)-1
%                 bordLoc = (rPos(:,i)<BC_x(j+1) & rPos(:,i)>=BC_x(j));
%                 SampChoice = pSamp<=BC(1,j);
%                 SampChoice2 = pSamp<=BC(2,j);
%                 newTrack(bordLoc & during & SampChoice & CW_States & nTurn(:,i)<3) = true;
%                 newTrack(bordLoc & during & SampChoice2 & CW_States & nTurn(:,i)>=3) = true;
%             end
%         end
        
        [state(newTrack,i)] = getNewState(state(newTrack,i-1),TP,before(newTrack),f(newTrack),df(newTrack));
        
        % if the fly is at the boundary, then it's in the boundary state
        state(bound & newTrack,i) = Bound_Ndx;
        newTrackState = state(:,i);
        % update the states that are not new tracks
        state(state(:,i)==0,i) = state(state(:,i)==0,i-1);
        
        % calculate time since inhibition
        durInhibition = (f < 0.01) & (abs(df) <1);%f<fBaseline;
        tt_inhibition(~durInhibition,i) = -1;
        if i>1
            tt_inhibition(durInhibition,i) = tt_inhibition(durInhibition,i-1)+1;% increase inhibition count by 1
        end
        ttInh = tt_inhibition(:,i);ttInh(~durInhibition) = nan;
        
        
        nTurn(newTrackState==ST_Ndx,i) = nTurn(newTrackState==ST_Ndx,i-1)+1;                % add to turns counter
    end
    
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % set up routines
    
    % cases 1-3 follow the new track routines
    cases = zeros(totFly,9);
    cases(:,1) = (newTrackState==1);
    cases(:,2) = (newTrackState==2);
    cases(:,3) = (newTrackState==3);
    
    % cases 4-5 follow the new boundary routines
    cases(:,4) = (newTrackState==4) & before;
    cases(:,5) = (newTrackState==4) & during;
    
    % cases 6-8 follow the continue movement routine
    cases(:,6) = ~newTrack & bound;
    cases(:,7) = ~newTrack & (state(:,i)==ST_Ndx); % stop and sharp turn have same continue movement routine
    cases(:,8) = ~newTrack & (state(:,i)==Stop_Ndx | state(:,i)==CW_Ndx);% curved walk has own routine
    
    % case 9 follow the move away from edge routine
    if i>1 
        if any(loc(:,i-1)==1)
            cases(:,9) = newTrack & bound & loc(:,i-1)==1;
        end
    end
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    % compute new routines (1-11)
    % new track routines
    for c = 1:3
        if any(cases(:,c))
            cc = cases(:,c)==true;
            
            f = nanmean(currSpk(:,max((i-hist-1),1):i-1),2);%nanmean(currSpk(:,max((i-6.*floor(ratio)),1):i),2);
            df = nanmean(dSpk(:,max((i-hist-1),1):i-1),2);%nanmean(dSpk(:,max((i-6.*floor(ratio)),1):i),2);
            df(isnan(df)) = 0;f(isnan(f)) = currSpk(1);% set nans to default 0 hz change and baseline spiking
            dfSmooth(:,i) = df;
            fSmooth(:,i) = f;
            df(abs(df)>150) = sign(df(abs(df)>150)).*150;
            
            durBaseline = cc & abs(df)<0.01 & abs(f-fBaseline)<0.01 & t_relFE(:,i)>0;% no change in frequency, frequency at baseline, and after first entry            
            tt_baseline(~durBaseline,i) = -1;
            if i>1
                tt_baseline(durBaseline,i) = tt_baseline(durBaseline,i-1)+1;% increase baseline count by 1
            end
            ttBase = tt_baseline(:,i);ttBase(~durBaseline) = -1;
            
            % sample kinematics and duration
            if c == ST_Ndx
                r = zeros(sum(cc),1);
                [r(:,1)] = sampleKinematics_raw040521(self,ST_spd,t_relFE(cc,i),f(cc),df(cc),ttBase(cc),ttInh(cc)./fs,1);
                [r(:,2)] = sampleKinematics_raw040521(self,ST_curv,t_relFE(cc,i),f(cc),df(cc),ttBase(cc),ttInh(cc)./fs,1);%1.5
                [r(:,3)] = sampleKinematics_raw040521(self,ST_dur,t_relFE(cc,i),f(cc),df(cc),ttBase(cc),nanVec(cc),1);
                r(:,3) = r(:,3);
                
                direction = sampleTurnBias040721(self,turnBias_ST,...
                    t_relFE(cc,i),f(cc),df(cc),ttBase(cc),...
                    xPos(cc,max(1,i-hist):i),yPos(cc,max(1,i-hist):i),...
                    ttInh(cc)./fs,r(:,2),r(:,1),r(:,3).*fs,'ST');%turnBias_ST
                
                currSpd(cc) = r(:,1);
                currCurv(cc) = r(:,2);%.*fs./originalFs
                oldDur(cc) = ceil(r(:,3).*fs);
                
            elseif c == CW_Ndx
                r = zeros(sum(cc),1);
                [r(:,1)] = sampleKinematics_raw040521(self,CW_spd,t_relFE(cc,i),f(cc),df(cc),ttBase(cc),ttInh(cc)./fs,1);
                [r(:,2)] = sampleKinematics_raw040521(self,CW_curv,t_relFE(cc,i),f(cc),df(cc),ttBase(cc),ttInh(cc)./fs,1);%1.5
                [r(:,3)] = sampleKinematics_raw040521(self,CW_dur,t_relFE(cc,i),f(cc),df(cc),ttBase(cc),nanVec(cc),1);
                
                direction = sampleTurnBias040721(self,turnBias_CW,...
                    t_relFE(cc,i),f(cc),df(cc),ttBase(cc),...
                    xPos(cc,max(1,i-hist):i),yPos(cc,max(1,i-hist):i),...
                    ttInh(cc)./fs,r(:,2),r(:,1)./fs./self.rBound./10,r(:,3).*fs,'CW');%turnBias_CW
                
                currSpd(cc) = r(:,1);
                currCurv(cc) = mod(r(:,2)./fs,360);%.*originalFs
                oldDur(cc) = ceil(r(:,3).*fs);
                
            elseif c == Stop_Ndx
                r = zeros(sum(cc),1);
                r(:,1) = sampleKinematics_raw040521(self,Stop_curv,t_relFE(cc,i),f(cc),df(cc),ttBase(cc),nanVec(cc),1);
                r(:,2) = sampleKinematics_raw040521(self,Stop_dur,t_relFE(cc,i),f(cc),df(cc),ttBase(cc),nanVec(cc),1);
                currSpd(cc) = 10.^-10;%need a small speed since curvature is normalized by speed (wrong when dividing by 0)
                oldDur(cc) = round(r(:,2).*fs);
                currCurv(cc) = r(:,1)./oldDur(cc);%mod(r(:,1),360);%.*fs./originalFs
                
                direction = sampleTurnBiasStops(self,turnBias_Stop,...
                    t_relFE(cc,i),f(cc),df(cc),ttBase(cc),...
                    xPos(cc,max(1,i-hist):i),yPos(cc,max(1,i-hist):i),currCurv(cc));
                
            end
            
            % correlated direction with last track
            if i>1
                currCurv(cc) = currCurv(cc).*direction;
            else
                currCurv(cc) = currCurv(cc).*sign(rand(sum(cc),1)-0.5);
            end
            
            currDur(cc) = oldDur(cc)-1;
            currSpd(cc) = currSpd(cc)./fs./self.rBound./10;% speed in normalized radial units/frame
            
            
            if (c == CW_Ndx) || (c == Stop_Ndx) % if curved walks, update all kinematics
                [spdAll(cc,i),curvAll(cc,i),phiAll(cc,i)] = updateGlobalKin(currSpd,currCurv,nanVec,cc);
            else % if stops or sharp turns, wait until mid point to update curvature
                [spdAll(cc,i),~,phiAll(cc,i)] = updateGlobalKin(currSpd,currCurv,nanVec,cc);
                curvAll(cc,i) = 0;
            end
            
            % update angle for curved walk only
            if c == CW_Ndx || c == Stop_Ndx
                currAng(cc) = updateCurrAng(allAng(:,i),curvAll(:,max(i-1,1)),curvAll(:,i),cc);% allAng(:,i) is the previous ang
            else
                currAng(cc) = allAng(cc,i);
            end
            
        end
    end
    allRout1 = any(cases(:,1:3),2);
    xPos(allRout1,i+1) = xPos(allRout1,i)+currSpd(allRout1).*cosd(currAng(allRout1));
    yPos(allRout1,i+1) = yPos(allRout1,i)+currSpd(allRout1).*sind(currAng(allRout1));

    % new boundary routines
    for c = 4:5
        if any(cases(:,c))
            cc = cases(:,c)==1;
            r = sampleKinematics_Boundary(self,c-3,sum(cc));%Joint
            
            currPhi(cc) = r(:,1);
            oldDur(cc) = round(r(:,2)).*ratio;%fs./originalFs;
            currDur(cc) = oldDur(cc)-1;
            [spdAll(cc,i),curvAll(cc,i),phiAll(cc,i)] = updateGlobalKin(nanVec,nanVec,currPhi,cc);
        end
    end
    % update the positional of the new boundary routine flies by rotating
    % their position by the angle theta
    allRout2 = find(cases(:,4) | cases(:,5));
    if ~isempty(allRout2)
        for j = 1:length(allRout2)
            k = allRout2(j);
            theta = currPhi(k);
            Arot = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
            xy = Arot * [xPos(k,i); yPos(k,i)];
            xy2 = Arot * xy;
            xPos(k,i+1:i+2) = [xy(1) xy2(1)];
            yPos(k,i+1:i+2) = [xy(2) xy2(1)];
        end
        
        currSpd(allRout2) = sqrt(diff(xPos(allRout2,i:i+1),[],2).^2+diff(yPos(allRout2,i:i+1),[],2).^2);
        currAng(allRout2) = myatan(diff(xPos(allRout2,i:i+1),[],2)',diff(yPos(allRout2,i:i+1),[],2)','degrees',2);
        tmpCurv = updateCurvatureFromTrack(xPos(:,max(i-2,1):i+2),yPos(:,max(i-2,1):i+2),allRout2);
        curvAll(allRout2,max(i-1,1):i) = tmpCurv(:,min(2,end-1):end-1);
        [spdAll(allRout2,i),~,~] = updateGlobalKin(currSpd,nanVec,currPhi,allRout2);
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % continue movement routines
    % boundary condition
    if any(cases(:,6))
        c = find(cases(:,6));
        for j = 1:length(c)
            k = c(j);
            theta = currPhi(k);
            % update the positional of the new boundary routine flies by 
            % rotating their position by the angle theta
            Arot = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
            xy = Arot * [xPos(k,i); yPos(k,i)];
            %xPos(k,i+1) = xy(1);
            %yPos(k,i+1) = xy(2);
            xy2 = Arot * xy;
            xPos(k,i+1:i+2) = [xy(1) xy2(1)];
            yPos(k,i+1:i+2) = [xy(2) xy2(1)];
        end
        currDur(c) = currDur(c)-1;
        
        currSpd(c) = sqrt(diff(xPos(c,i:i+1),[],2).^2+diff(yPos(c,i:i+1),[],2).^2);
        currAng(c) = myatan(diff(xPos(c,i:i+1),[],2)',diff(yPos(c,i:i+1),[],2)','degrees',2);
        tmpCurv = updateCurvatureFromTrack(xPos(:,max(i-2,1):i+2),yPos(:,max(i-2,1):i+2),c);
        curvAll(c,max(i-1,1):i) = tmpCurv(:,min(2,end-2):end-1);
        [spdAll(c,i),~,phiAll(c,i)] = updateGlobalKin(currSpd,nanVec,currPhi,c);
    end
    % sharp turn/stop only updates angle halfway through the duration
    if any(cases(:,7))
        c = cases(:,7)==1;
        halfDur = (floor(oldDur./2)+1);
        %halfDur = max(oldDur-1,1);% at the start
       
        [spdAll(c,i),~,phiAll(c,i)] = updateGlobalKin(currSpd,currCurv,nanVec,c);
        curvAll(c,i) = 0;
        [~,curvAll(c & currDur==halfDur,i),~] = updateGlobalKin(currSpd,currCurv,nanVec, c & currDur==halfDur);
        currAng(c) = updateCurrAng(allAng(:,i),curvAll(:,max(i-1,1)),curvAll(:,i),c);
        
        xPos(c,i+1) = xPos(c,i)+currSpd(c).*cosd(currAng(c));
        yPos(c,i+1) = yPos(c,i)+currSpd(c).*sind(currAng(c));
        currDur(c) = currDur(c)-1;
    end
    % curved walk continuously update movement position and angle
    if any(cases(:,8))
        c = cases(:,8)==1;
        
        [spdAll(c,i),curvAll(c,i),phiAll(c,i)] = updateGlobalKin(currSpd,currCurv,nanVec,c);
        
        %currAng(c) = currAng(c)+currCurv(c);
        currAng(c) = updateCurrAng(allAng(:,i),curvAll(:,max(i-1,1)),curvAll(:,i),c);
        xPos(c,i+1) = xPos(c,i)+currSpd(c).*cosd(currAng(c));
        yPos(c,i+1) = yPos(c,i)+currSpd(c).*sind(currAng(c));
        currDur(c) = currDur(c)-1;
        
        
    end
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    % move away from edge routine
    if any(cases(:,9))
        c = cases(:,9)==1;
        %move the fly away from the boundary
        xPos(c,i+1) = xPos(c,i)./1.02;
        yPos(c,i+1) = yPos(c,i)./1.02;
        xPos(c,i+2) = xPos(c,i+1)./1.02;
        yPos(c,i+2) = yPos(c,i+1)./1.02;
        
        
        % reorient to face directly inwards 
        currAng(c) = myatan(xPos(c,i+1)',yPos(c,i+1)','degrees',2)+180;%+20*(rand(1)-0.5);
        currDur(c) = 0;
        oldDur(c) = 0;
        turn(c) = false;
        
        currSpd(c) = sqrt(diff(xPos(c,i:i+1),[],2).^2+diff(yPos(c,i:i+1),[],2).^2);
        tmpCurv = updateCurvatureFromTrack(xPos(:,i-2:i+2),yPos(:,i-2:i+2),find(c));
        curvAll(c,i-1:i) = tmpCurv(:,2:end-1);
        [spdAll(c,i),~,phiAll(c,i)] = updateGlobalKin(currSpd,nanVec,nanVec,c);
    end
    %----------------------------------------------------------------------
    
    allAng(:,i+1) = currAng;
    % update x,y, and r position for points that move out of the arena
    % bounds
    xTemp = xPos(:,i+1);
    yTemp = yPos(:,i+1);
    
    rTemp = sqrt(xPos(:,i+1).^2+yPos(:,i+1).^2);
    rTemp(rTemp<1) = 1;
    
    xPos(:,i+1) = xTemp./rTemp;
    yPos(:,i+1) = yTemp./rTemp;
    rPos(:,i+1) = sqrt(xPos(:,i+1).^2+yPos(:,i+1).^2);
    
    newTrack(:) = false;
    progressbar(i/(len-1))
    
    if mod(i,1*60*fs)==1
        disp([num2str(i) '/' num2str(len)])
    end
end
curvAll = curvAll.*pi./180;% convert curvature to radians since emp is in radians
xPos = xPos(:,1:len);
yPos = yPos(:,1:len);
rPos = rPos(:,1:len);

% calculate first entry
firstEntry = ones(totFly,1)*len;
for i = 1:totFly
    tmpFE = find(rPos(i,ceil(lightOnTime):len)<border,1);
    if ~isempty(tmpFE)
        firstEntry(i) = tmpFE+floor(lightOnTime)-1;
    end
end

% plot the first 6 tracks
figure;set(gcf,'Position',[1 40 585 600])
for i = 1:min(6,totFly)
    subplot(3,2,i);
    plot(xPos(i,1:firstEntry(i)),yPos(i,1:firstEntry(i)),'g');hold on
    plot(xPos(i,firstEntry(i):end),yPos(i,firstEntry(i):end),'r');
    plot(xPos(i,1),yPos(i,1),'k*')
    axis([-1 1 -1 1])
end

% assign output structure
synthFlys.x = xPos.*4;
synthFlys.y = yPos.*4;
synthFlys.r = rPos.*4;
synthFlys.firstEntry = firstEntry;
synthFlys.lightOn = lightOnTime.*ones(totFly,1);
synthFlys.spk = currSpk;
end

% update speed, curv, phi matrix
function [spdAll,curvAll,phiAll] = updateGlobalKin(currSpd,currCurv,currPhi,c)
spdAll = currSpd(c);
curvAll = currCurv(c);
phiAll = currPhi(c);
end


% update speed, curv, phi matrix
function [currAng] = updateCurrAng(prevAng,prevCurv,currCurv,c)
currAng = prevAng(c)+(prevCurv(c)+currCurv(c))./2;
end

function [currCurv] = updateCurvatureFromTrack(xPos,yPos,c)
currCurv = nan(length(c),size(xPos,2)-1);
for j = 1:length(c)
    k = c(j);
    Vertices = horzcat(xPos(k,:)',yPos(k,:)');
    %get normal vectors
    N=LineNormals2D(Vertices);
    
    theta=acos(N(:,1));
    opp = (N(:,1)>0 & N(:,2)<0) | (N(:,1)<0 & N(:,2)<0);
    theta(opp) = -theta(opp);
    theta = unwrap(theta);
    curv = diff(theta).*180./pi;
    currCurv(j,:) = curv;
end

end
