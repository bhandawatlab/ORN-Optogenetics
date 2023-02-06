function [self] = getBaselineKinematicsEarlyLate(self,border,condNdx,key)
% getBaselineKinematicsEarlyLate  calculates the kinematics when flies
%   are at baseline firing. This is split into before first entry, the
%   first 2 state instances after each exit, and lateral state instances
%   after each exit
%
%   Inputs: self = fly object
%           border = radial position of light border (in cm)
%           condNdx = fly x time matric of the general firing rate regime
%           the fly is in (e.g. 1=before, 2=below, 3=baseline, and 4=above)
%           key = cell array key for condNdx
%
%   Output: self = fly object updated with the baseline kinematics
%   

stateDistEarly = cell(4,1);stateDistLate = cell(4,1);stateDistBefore = cell(4,1);

baseline = getBaselineKinematics2(self,border,condNdx,key);

allTracks = baseline.during;
allStates = allTracks.state;

earlyNdx = 1:2;
lateNdx = earlyNdx(end)+1:10;%size(allStates,2);

nStates = numel(self.states.key);

for state = 1:nStates
    stateDistEarly{state}.curv = [];
    stateDistEarly{state}.spd = [];
    stateDistEarly{state}.dur = [];
    
    stateDistLate{state}.curv = [];
    stateDistLate{state}.spd = [];
    stateDistLate{state}.dur = [];
    
    stateDistBefore{state}.curv = [];
    stateDistBefore{state}.spd = [];
    stateDistBefore{state}.dur = [];
end

% [stateDistBefore] = getKinematics(self,baseline.before.state,baseline.before,stateDistBefore,[1:1000]);
% % sanity check
% if ~isempty(self.model.params)
%     % may contain 1 extra track in some flies (<1% difference)
%     tmp2 = [nanmean(stateDistBefore{1,1}.spd), nanmean(self.model.params{1,1}.baselineDat');...
%         nanmean(stateDistBefore{2,1}.spd), nanmean(self.model.params{1,2}.baselineDat');...
%         nanmean(stateDistBefore{1,1}.curv), nanmean(self.model.params{1,3}.baselineDat');...
%         nanmean(stateDistBefore{2,1}.curv), nanmean(self.model.params{1,4}.baselineDat');...
%         nanmean(stateDistBefore{1,1}.dur), nanmean(self.model.params{1,5}.baselineDat');...
%         nanmean(stateDistBefore{2,1}.dur), nanmean(self.model.params{1,6}.baselineDat');...
%         nanmean(stateDistBefore{3,1}.curv), nanmean(self.model.params{1,7}.baselineDat');...
%         nanmean(stateDistBefore{3,1}.dur), nanmean(self.model.params{1,8}.baselineDat')];
% end
[stateDistEarly] = getKinematics(self,allStates,allTracks,stateDistEarly,earlyNdx);
[stateDistLate] = getKinematics(self,allStates,allTracks,stateDistLate,lateNdx);


% fit raw values to distribution and save to the current fly object if the
% model model parameters have already been initialized
if ~isempty(self.model)
    for i = 1:numel(self.model.params)
        kin = self.model.params{i}.state.kin;
        state = self.model.params{i}.state.state;
        
        if strcmpi(kin,'totCurv') || strcmpi(kin,'avgCurv')
            kin = 'curv';
        end
        
        l = cellfun(@(c)strcmp(c,self.states.key),{state},'UniformOutput',false);
        ndx = find(l{1});
        
        % early baseline after first entry
        x = stateDistEarly{ndx}.(kin);
        if ~isempty(x)
            parmhat =  lognfit(x+eps);
        else
            % set to before first entry if no data
            parmhat =  self.model.params{i}.baseline;
        end
        self.model.params{i}.baselineEarly = parmhat;
        
        % later baseline after first entry
        x = stateDistLate{ndx}.(kin);
        if ~isempty(x)
            parmhat =  lognfit(x+eps);
        else
            % set to early baseline if no data
            parmhat =  self.model.params{i}.baselineEarly;
        end
        self.model.params{i}.baselineLate = parmhat;
    end
end


end


function [stateDist] = getKinematics(self,allStates,allTracks,stateDist,ndx)

nStates = numel(self.states.key);
for fly = 1:self.nFly
    currFlyStates = squeeze(allStates(fly,:,:));
    goodTracks = find(~all(isnan(currFlyStates),2));
    
    for track = 1:numel(goodTracks)
        currTrackStates = currFlyStates(goodTracks(track),:);
        currTrackStates(isnan(currTrackStates)) = [];
        
        for state = 1:nStates
            [startNdx,endNdx,type] = startEndSeq(currTrackStates==state);
            startNdx = startNdx(type);endNdx = endNdx(type);
            len = endNdx-startNdx+1;
            startNdx(len<3) = [];endNdx(len<3) = [];
            
            for i = 1:numel(startNdx)
                if strcmpi(self.states.key{state},'sharp turns')
                    tmpCurve = abs(nansum(allTracks.curv(fly,goodTracks(track),startNdx(i):endNdx(i))));
                    tmpCurve = mod(tmpCurve,360);
                elseif strcmpi(self.states.key{state},'curved walks')
                    tmpCurve = nanmean((allTracks.curv(fly,goodTracks(track),startNdx(i):endNdx(i)))).*self.fs;
                elseif strcmpi(self.states.key{state},'stops')
                    tmpCurve = abs(nansum(allTracks.curv(fly,goodTracks(track),startNdx(i):endNdx(i))));
                    tmpCurve = mod(tmpCurve,360);
                elseif strcmpi(self.states.key{state},'boundary')
                    %tmpCurve = diff(abs(allTracks.curv(fly,goodTracks(track),[startNdx(i),endNdx(i)])));
                    tmpCurve = nansum(allTracks.curv(fly,goodTracks(track),startNdx(i):endNdx(i)));
                end
                tmpSpd = nanmean(allTracks.spd(fly,goodTracks(track),startNdx(i):endNdx(i)));
                tmpDur = (endNdx(i)-startNdx(i)+1)./self.fs;
                
                if any(i == ndx)
                    stateDist{state}.curv = [stateDist{state}.curv, abs(tmpCurve)];
                    stateDist{state}.spd = [stateDist{state}.spd, tmpSpd];
                    stateDist{state}.dur = [stateDist{state}.dur, tmpDur];
                end
            end
        end
    end
end

end