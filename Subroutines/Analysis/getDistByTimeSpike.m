function [self] = getDistByTimeSpike(self,condNdx,key,opts,plotFig)
% getDistByTimeSpike  wrapper function for calculating the transition 
%   probability, turn optimality, and kinematics in the KNN space
%
%   Inputs: self = fly object
%           condNdx = fly x time matric of the general firing rate regime
%           the fly is in (e.g. 1=before, 2=below, 3=baseline, and 4=above)
%           key = cell array key for condNdx
%           opts = meta information about how to set up KNN space
%           plotFig = true/false on whether to plot figures
%
%   Output: self = fly object updated with the KNN space
%   

%% Set up the data
% lump all kinematic changes caused by ORN activity into one category
spkTrigKinNdx = find([cellfun(@(x) strcmp(x,'below'),key,'UniformOutput',true) |...
    cellfun(@(x) strcmp(x,'above'),key,'UniformOutput',true)]);
% update the condition matrix into before, during baseline, and during
% nonbaseline
condNdx2 = zeros(size(condNdx));
for i = 1:numel(spkTrigKinNdx)
    condNdx2(condNdx==spkTrigKinNdx(i)) = 1;
end
% get all kinematics and other values separated by locomotor state and
% sensory condition (i.e 1 = before, 2 = non-baseline, 3 = during baseline)
[s, keyNdx] = getKinematicConditions(self,condNdx2,opts);

% set x (delta firing rate),y (firing rate),and z (time) points to get KNN estimation 
xGrid = opts.xGrid;
yGrid = opts.yGrid;
zGrid = opts.zGrid;%[0:5:180].*f_orco.fs;
% set the radii of the bounding ellipse for KNN estimations
ratio = opts.ratio;%[30, 10, 20.*f_orco.fs];

% Define each state kinematic
allCond{1}.state = 'sharp turns';allCond{1}.kin= 'spd';
allCond{2}.state = 'curved walks';allCond{2}.kin= 'spd';
allCond{3}.state = 'sharp turns';allCond{3}.kin= 'totCurv';
allCond{4}.state = 'curved walks';allCond{4}.kin= 'avgCurv';
allCond{5}.state = 'sharp turns';allCond{5}.kin= 'dur';
allCond{6}.state = 'curved walks';allCond{6}.kin= 'dur';
allCond{7}.state = 'stops';allCond{7}.kin= 'totCurv';
allCond{8}.state = 'stops';allCond{8}.kin= 'dur';


%% get the probability of transitioning between states
fprintf('Calculating Transition Probability\n');
% set the max distance and number of data points for the KNN estimation
if isempty(opts.Thresh)
    Thresh = 1.5;
else
    Thresh = opts.Thresh;
end
if isempty(opts.K)
    K = 128; % note that this K is larger because we have more data (not spread out across time)
else
    K = opts.K;
end
for stateN = 1:numel(self.states.key)
    state = self.states.key(stateN);
    try
        [newStateProb,XX,YY,newState,newStateProbBaseline,newStateProbBaselineDuring]...
            = getKNNTransProb(self, s, keyNdx, xGrid, yGrid, ...
            ratio, Thresh, K, state, plotFig);
        self.model.TP.before(stateN,1:numel(newStateProbBaseline)) = newStateProbBaseline;
        self.model.TP.during(stateN,1:numel(newStateProb)) = newStateProb;
        self.model.TP.duringBaseline(stateN,1:numel(newStateProbBaselineDuring)) = newStateProbBaselineDuring;
        self.model.TP.key = self.states.key;
        self.model.TP.XX = XX;
        self.model.TP.YY = YY;
    catch
        self.model.TP.before(stateN,1:(numel(self.states.key)-1)) = nan(1,(numel(self.states.key)-1));
        self.model.TP.during(stateN,1:(numel(self.states.key)-1)) = [];
        self.model.TP.key = self.states.key;
%         self.model.TP.XX = XX;
%         self.model.TP.YY = YY;
    end
end

%% get turn optimality for sharp turns, curved walks, and stops
fprintf('Calculating Turn Optimality\n');
% set the max distance and number of data points for the KNN estimation
if isempty(opts.Thresh)
    Thresh = 1.5;
else
    Thresh = opts.Thresh;
end
if isempty(opts.K)
    K = 64;
else
    K = opts.K;
end
% calculate turn bias only on states 1-3 (ST, CW, and stops, not boundary)
stateTB = [1:3];state = self.states.key(stateTB);
% calculate turn optimality
[turnBias,turnBiasRaw,XX,YY,turnBias_baseline,turnBias_baseline_during,turnBias_during_baseline_byRPos] = ...
    getKNNProbTurnIn(self, s, keyNdx, xGrid, yGrid, ...
    zGrid,ratio, Thresh, K, state, opts.border, self.id, plotFig);
% update turn optimality fields in the fly object
self.model.TurnBias.before = turnBias_baseline;
self.model.TurnBias.during = turnBias;
self.model.TurnBias.duringRaw = turnBiasRaw;
self.model.TurnBias.during_baseline = turnBias_baseline_during;
self.model.TurnBias.during_baseline_byRPos = turnBias_during_baseline_byRPos;
self.model.TurnBias.key = self.states.key(stateTB);
% update the KNN parameters used to calculate the turn optimality
self.model.TurnBias.XX = XX{1}(:,:,1);
self.model.TurnBias.YY = YY{1}(:,:,1);
self.model.TurnBias.tt = zGrid./self.fs;
self.model.TurnBias.dt = ratio(3)./self.fs;

%% inhibition turn bias and turn optimality
% for cond = 3:4
%     Thresh = 1.5; K = 64;
%     getKNNKin(f_orco, s, keyNdx, xGrid, yGrid, zGrid, ratio, ...
%         Thresh, K, allCond{cond}.state, allCond{cond}.kin, f_orco.id, true);%plotFig
% end

%% optimize for K and window size for kinematics
% currFig = get(gcf,'Number');
% figure(currFig+1);set(gcf,'Position',[2 42 838 924])
% figure(currFig+2);set(gcf,'Position',[2 42 838 924])
% fprintf('Optimizing for K and Thresh\n');
% for cond = 1:numel(allCond) 
%     % for plotting purposes, these were the K-values used
%     if cond <7
%         Thresh = 1.5; K = 64;
%     else
%         Thresh = 1; K = 64;
%     end
%     % optimize for K and window size for kinematics
%     estimationSuccess = getParams(self, s, keyNdx, xGrid, yGrid, zGrid, ratio, allCond{cond}.state,...
%         allCond{cond}.kin, Thresh, K, opts.border, currFig+[1,2],[4,2,cond]);
%     if estimationSuccess
%         figure(currFig+1);
%         title([allCond{cond}.state ' ' allCond{cond}.kin])
%         figure(currFig+2);
%         title([allCond{cond}.state ' ' allCond{cond}.kin])
%     end
% end

%% KNN estimation of kinematics
for cond = 1:numel(allCond)
    fprintf('Generating KNN estimations, %d / %d\n',[cond, numel(allCond)]);
    % use some previously found K/window size value since the optimization 
    % takes time
    if isempty(opts.Thresh)
        if cond <7
            Thresh = 1.5;
        else
            Thresh = 1;
        end
    else
        Thresh = opts.Thresh;
    end
    if isempty(opts.K)
        if cond <7
            K = 64;
        else
            K = 64;
        end
    else
        K = opts.K;
    end
    %----------------------------------------------------------------------
    
    % Get the empirical distribution (logN) parameters as a function of a
    % grid of df and f and time since first entry
    try
    [fitFun,pHat,XX,YY,pHatBase,CorrMat1,CorrMat2,qBounds1,qBounds2,...
        baselineDat,allDat,allSpk,alldSpk,allTimeFE,allTime,allFly,...
        allFlyBaseline,turnBias,nLL,AIC,BIC,distTypes,maskedDat]...
        = getKNNKin(self, s, keyNdx, xGrid, yGrid, zGrid, ratio, ...
        Thresh, K, allCond{cond}.state, allCond{cond}.kin, opts.border, self.id, true, false);
    catch
        disp('failed to estimate KNN space')
    end
        
    
    % reshape the mu and sigma of the logN into dF x F x time matrix
    val = reshape(pHat(:,1),size(XX));
    val2 = reshape(pHat(:,2),size(XX));% std
    nanSlices = (squeeze(sum(~isnan(val),[1,2]))<20);
    
    %if all slices don'e have enought data points, then just set to
    %baseline
    if all(nanSlices)
        for slice = 1:numel(nanSlices)
            val(:,:,slice) = pHatBase(1);
            val2(:,:,slice) = pHatBase(2);
        end
    % if any time slice doesn't have enough data points
    elseif sum(nanSlices)>0
        % find the slices before and after that has enough data points
        [closestPtBef,closestPtAft,idxB] = findBeforeAfter(find(nanSlices),find(~nanSlices),'both');
        allNanSlice = find(nanSlices);
        
        closestPtBef(isnan(closestPtBef)) = 1;
        closestPtAft(isnan(closestPtAft)) = numel(nanSlices);
        % average the mu and sigma values based on the closest slices
        % before and after
        for slice = 1:sum(nanSlices)
            val(:,:,allNanSlice(slice)) = nanmean(val(:,:,[closestPtBef(slice),closestPtAft(slice)]),3);
            val2(:,:,allNanSlice(slice)) = nanmean(val2(:,:,[closestPtBef(slice),closestPtAft(slice)]),3);
        end
    end
    
    % for each time slice
    for slice = 1:numel(zGrid)
        % interpret the surface by smoothing and filling in any part of the
        % space that do not have enough data to estimate an KNN
        % distribution as the closest vertix that does.
        if ~all(isnan(val(:,:,slice)),'all')
            [g2,~,~,~] = interpretLandscape(val(:,:,slice)',XX(:,:,slice)',YY(:,:,slice)',2,false);
            [g2_std,~,~,~] = interpretLandscape(val2(:,:,slice)',XX(:,:,slice)',YY(:,:,slice)',2,false);

            % add model parameters to the fly object (kinematics)
            self.model.params{cond}.g2{1,slice} = g2;
            self.model.params{cond}.g2_std{1,slice} = g2_std;
        else
            self.model.params{cond}.g2{1,slice} = nan;
            self.model.params{cond}.g2_std{1,slice} = nan;
        end
    end
    % add model parameters to the fly object (kinematics)
    self.model.params{cond}.val = val;
    self.model.params{cond}.val2 = val2;
    self.model.params{cond}.state =  allCond{cond};
    self.model.params{cond}.fitFun = fitFun;
    self.model.params{cond}.baseline = pHatBase;
    % add model parameters to the fly object (all empirical data)
    self.model.params{cond}.baselineDat = baselineDat;
    self.model.params{cond}.allSpk = allSpk;
    self.model.params{cond}.alldSpk = alldSpk;
    self.model.params{cond}.allDat = allDat;
    self.model.params{cond}.allTime = allTime;
    self.model.params{cond}.allTimeFE = allTimeFE;
    self.model.params{cond}.allFly = allFly;
    self.model.params{cond}.allFlyBaseline = allFlyBaseline;
    self.model.params{cond}.turnBias = turnBias;
    self.model.params{cond}.KNN_rawData = maskedDat;
    
    % add model parameters to the fly object (KNN parameters)
    self.model.params{cond}.KNN.Thresh = Thresh;
    self.model.params{cond}.KNN.K = K;
    self.model.params{cond}.KNN.ratio = ratio;
    self.model.params{cond}.KNN.Time = zGrid;
    self.model.params{cond}.KNN.tt = zGrid./self.fs;
    
    % add model fits
    self.model.params{cond}.KNN.nLL = nLL;
    self.model.params{cond}.KNN.AIC = AIC;
    self.model.params{cond}.KNN.BIC = BIC;
    self.model.params{cond}.KNN.DistTypes = distTypes;
    
    % add model parameters to the fly object (correlations and bounds)
    self.model.params{cond}.corrBaseline = CorrMat1;
    self.model.params{cond}.corr = CorrMat2;
    self.model.params{cond}.qBoundsBaseline = qBounds1;
    self.model.params{cond}.qBounds = qBounds2;
    if strcmpi(allCond{cond}.state,'stops')
        self.model.params{cond}.corrRow = [{'curv'};{'dur'}];
    else
        self.model.params{cond}.corrRow = [{'spd'};{'curv'};{'dur'}];
    end
end
fprintf('Done generating KNN distributions\n')

end

function [s, keyNdx] = getKinematicConditions(f_orco,condNdx,opts)
    condNdx3 = zeros(size(f_orco.spk));
    condNdx3(condNdx<6) = 1;
    condNdx3(condNdx==0) = 3;% not an condition case
    contIdx = 2;% index of the control group (before light on period)
    fe = f_orco.getFirstEntry('H',opts.border);
    for i = 1:f_orco.nFly
        condNdx3(i,1:fe(i)-1) = contIdx;% before case
    end
    key = cellstr(num2str(unique(condNdx3)));thresh = 0;%0.1
    %[s, keyNdx] = f_orco.getKinematicsCond(condNdx3,key,thresh);
    [s, keyNdx] = f_orco.getKinematicsCond(condNdx3,key,thresh,opts.timeInterval);
end











