function [estimationSuccess] = getParams(f_orco, s, keyNdx, xGrid, yGrid, zGrid, ratio, ...
    state, kin, Thresh, K, border, figN,sbpltN)
% getParams  does a grid search to determine model params. Note that this
%   does not automaticaly look for the elbow point. The user should look at
%   the cost function to make a decision.
%
%   Inputs: f_orco = fly object
%           s = structure of descriptions for each state (e.g. average 
%               curvature, average speed, time ndx, prior ORN firing rate, 
%               etc). Each field (description) is a fly x condition cell
%               array. Each cell is a m x 1 vector of descriptors for each
%               state instance
%           keyNdx = Cell array describing the condition (each column) of
%               the cell arrays 
%           xGrid = x grid locations of KNN estimations (delta firing rate)
%           yGrid = y grid locations of KNN estimations (firing rate)
%           zGrid = z grid locations of KNN estimations (time since first entry)
%           ratio = the radii of the bounding ellipse for KNN estimations
%           state = index of the state that you are looking at
%           kin = string stating what kinematic feature you are looking at
%               (i.e. 'spd','avgCurv','totCurv', or 'dur')
%           Thresh = plotting purposes only (proposed thresh) - not
%               relevent for grid search
%           K = plotting purposes only (proposed K) - not
%               relevent for grid search
%           border = radial light border (in cm)
%           figN = figure number to plot results to
%           sbpltN = subplot number to plot results to ([row, column, index])
%
%   Output: estimationSuccess = true/false whether the analysis was
%           successful
%   

% calculate time from first entry for the beginning of each state
fe = f_orco.getFirstEntry('H',border);
for fly = 1:size(s.ndx,1)
    for j = 1:size(s.ndx,2)
        s.timeFromFe{fly,j} = s.ndx{fly,j}-fe(fly);
    end
end

% index in the s structure for the current locomotor state and after first
% entry (non-baseline)
l = cellfun(@(c)strcmp(c,keyNdx),{state},'UniformOutput',false);
l2 = cellfun(@(c)strcmp(c,keyNdx),{'1'},'UniformOutput',false);
typeIdx = any(l{1},2) & any(l2{1},2);

% index based on the kinematics 
if strcmpi(kin, 'spd')
    allDat = abs(cell2mat(s.avgSpd(:,typeIdx)));
elseif strcmpi(kin, 'avgCurv')
    allDat = abs(cell2mat(s.avgCurv(:,typeIdx))*180./pi).*f_orco.fs;%deg/s
elseif strcmpi(kin, 'totCurv')
    allDat = abs(cell2mat(s.totCurv(:,typeIdx))*180./pi);
    allDat = mod(allDat,360);
    badNdx = allDat<0.01;allDat(badNdx) = 0.01;
elseif strcmpi(kin, 'dur')
    allDat = abs(cell2mat(s.dur(:,typeIdx)))./f_orco.fs;% in seconds
end

% get the firing rate, delta firing rate, and time since first entry
% associated with each kinematic data point
allSpk = cell2mat(s.spkStart(:,typeIdx));
alldSpk = cell2mat(s.dSpkStart(:,typeIdx));
allTime = max(0,cell2mat(s.timeFromFe(:,typeIdx)));


%--------------------------------------------------------------------------
% set up the meshgrid
X = [alldSpk, allSpk, allTime]./ratio;
[XX,YY,ZZ] = ndgrid(xGrid,yGrid,zGrid);
Y = [XX(:),YY(:),ZZ(:)]./ratio;

% list of K and thresholds to grid search through
K2Cons = 2.^(4:1:8);
thresh2Cons = 0.5:0.5:5;
nPtsMin = 15;%minimum number of points to consider

% make sure that K doesn't exceed the maximum number of data points in the
% dataset
nPtMax = size(X,1);
K2Cons(K2Cons>nPtMax) = nPtMax;

% perform KNN search with the maximum K
[idx, D] = knnsearch(X,Y,'K',max(K2Cons));
maskedDat = log(allDat(idx));% lognormal
% look through and compute cost function
cost = zeros(numel(K2Cons),numel(thresh2Cons));
for k = 1:numel(K2Cons)
    for t = 1:numel(thresh2Cons)
        % set any data points outside of threshold to nans
        mask = D(:,1:K2Cons(k))<thresh2Cons(t);
        maskedDatTmp = maskedDat(:,1:K2Cons(k));
        maskedDatTmp(~mask) = nan;
        % calculate the standard dev. and SEM of each 
        stdVals = nanstd(maskedDatTmp,[],2);
        semVals = stdVals./sqrt(sum(mask,2));
        semVals(sum(mask,2)<nPtsMin) = nan;
        cost(k,t) = nanmean(semVals)./K2Cons(k)./thresh2Cons(t);
    end
end
if ~all(isnan(cost),'all')
    % plot the cost function
    figure(figN(1));subplot(sbpltN(1),sbpltN(2),sbpltN(3));
    [kk,tt] = ndgrid(K2Cons,thresh2Cons);
    surf(log2(kk),tt,cost);hold on;scatter3(log2(K),Thresh,cost(kk==K & tt==Thresh),1000,'.r')
    xlabel('K');ylabel('Threshold')
    zlabel('cost function (mean(sem)/(K*thresh))');view(135,45)
    %zlabel('cost function (mean(sem))');view(135,45)
    colormap(parula)
    yticks(thresh2Cons(2:2:end));ylim([0.5 5]);xlim([4 8])
    set (gca, 'XTickLabel', strcat('2^{',num2str(log2(K2Cons(:))),'}'));

    figure(figN(2));subplot(sbpltN(1),sbpltN(2),sbpltN(3));
    imagesc(thresh2Cons,log2(K2Cons),cost);hold on;scatter(Thresh,log2(K),1000,'.r')
    xlabel('Threshold');ylabel('K');
    colormap(parula);colorbar;
    xticks(thresh2Cons(2:2:end));xlim([0.5 5]);ylim([3.5 8.5])
    set (gca, 'YTickLabel', strcat('2^{',num2str(log2(K2Cons(:))),'}'));
    estimationSuccess = true;
else
    disp('Not enough data points to estimate K and Threshold for KNN space ')
    estimationSuccess = false;
end
end











