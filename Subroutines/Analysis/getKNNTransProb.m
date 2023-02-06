function [newStateProb,XX,YY,newState,newStateProbBaseline,newStateProbBaselineDuring]...
    =  getKNNTransProb(f_orco,s,keyNdx,xGrid,yGrid,ratio,Thresh,K,state,plotFig)
% getKNNTransProb  calculates the transition probability in the KNN space
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
%           ratio = the radii of the bounding ellipse for KNN estimations
%           Thresh = max distance for the KNN estimations
%           K = number of data points for the KNN estimations
%           state = index of the state that you are transitioning out of
%           plotFig = true/false of whether to plot analysis figure
%
%   Output: newStateProb = cell array of transition probability in the KNN
%               space for input state type. Note that each cell is a
%               function that takes in dF and F coordinates.
%           XX = dF coordinates in KNN space
%           YY =  F coordinates in KNN space
%           newState = same as newStateProb, but the non-interpolated
%               version, and is saved as a matrix rather than a function.
%           newStateProbBaseline = vector of transition probabilities
%               before first entry for the input state type
%           newStateProbBaselineDuring = vector of transition probabilities
%               after first entry, but when the fly is experiencing 
%               baseline firing rate for the input state type
%   

l = cellfun(@(c)strcmp(c,keyNdx),{state},'UniformOutput',false);
l2 = cellfun(@(c)strcmp(c,keyNdx),{'1'},'UniformOutput',false);%during non-baseline
l3 = cellfun(@(c)strcmp(c,keyNdx),{'2'},'UniformOutput',false);%before
l4 = cellfun(@(c)strcmp(c,keyNdx),{'3'},'UniformOutput',false);%during baseline
typeIdx = any(l{1},2) & any(l2{1},2);
baselineIdx = any(l{1},2) & any(l3{1},2);
during_baselineIdx = any(l{1},2) & any(l4{1},2);

newState_before = cell2mat(s.nextState(:,baselineIdx));
newState_during = cell2mat(s.nextState(:,typeIdx));
newState_duringBaseline = cell2mat(s.nextState(:,during_baselineIdx));
allSpk = cell2mat(s.spkEnd(:,typeIdx));
alldSpk = cell2mat(s.dSpkEnd(:,typeIdx));


% remove boundary
newState_before = newState_before(newState_before<4);
allSpk = allSpk(newState_during<4);
alldSpk = alldSpk(newState_during<4);
newState_during = newState_during(newState_during<4);
newState_duringBaseline = newState_duringBaseline(newState_duringBaseline<4);

% remove tracks without a previous track (initial tracks)
newState_before(isnan(newState_before)) = [];
allSpk(isnan(newState_during)) = [];
alldSpk(isnan(newState_during)) = [];
newState_during(isnan(newState_during)) = [];
newState_duringBaseline(isnan(newState_duringBaseline)) = [];

%get before transition probability
transitions = unique(newState_before);
newStateProbBaseline = zeros(1,numel(transitions));
for n = 1:numel(transitions)
    newStateProbBaseline(transitions(n)) = ...
        sum(newState_before==transitions(n))./numel(newState_before);
end

% get during baseline transition probability
transitions = unique(newState_duringBaseline);
newStateProbBaselineDuring = zeros(1,numel(transitions));
for n = 1:numel(transitions)
    newStateProbBaselineDuring(transitions(n)) = ...
        sum(newState_duringBaseline==transitions(n))./numel(newState_duringBaseline);
end

% get during KNN transition probability
[newStateProb,XX,YY,newState] = geTransitionProb(newState_during,alldSpk,allSpk,xGrid,yGrid,ratio,Thresh,K,newStateProbBaseline);

transitions = unique(newState_during);
if plotFig
    figure;set(gcf,'Position',[1 40 290 600])
    for n = 1:numel(transitions)
        subplot(numel(transitions),1,n);
        surf(XX,YY,newStateProb{transitions(n)}(XX,YY));view(45,45)%view(2);
        xlim([min(xGrid) max(xGrid)]);ylim([min(yGrid) max(yGrid)]);zlim([0 1])
        title(['P of transition into ' f_orco.states.key{transitions(n)}])
        
        suptitle(state)
    end
end

end

function [Vq,XX,YY,newStateAll] = geTransitionProb(newState,alldSpk, allSpk,xGrid,yGrid,ratio,Thresh,K,newStateProbBaseline)

X2 = [alldSpk, allSpk]./ratio(1:2);
[XX,YY] = meshgrid(xGrid,yGrid);
Y2 = [XX(:),YY(:)]./ratio(1:2);
[idx2, D2] = knnsearch(X2,Y2,'K',K);

newStateDat = newState(idx2);
mask2 = D2<Thresh;
newStateDat(~mask2) = nan;
newStateDat(sum(mask2,2)<30,:) = nan;
mask2(sum(mask2,2)<30,:) = 0;


if sum(sum(mask2,2)~=0)>(0.05*size(mask2,1))
    transitions = unique(newState);
    for n = 1:numel(transitions)
        pStop = sum(newStateDat==transitions(n),2)./sum(mask2,2);
        vals = reshape(pStop,size(XX));
        
        % set any nan values to the closest non-nan value
        [noDatr, noDatc] = find(isnan(vals));
        [datr, datc] = find(~isnan(vals));
        D = pdist2([noDatr, noDatc],[datr, datc]);
        [~,nearestPtNdx]=min(D,[],2);
        for i = 1:numel(noDatr)
            vals(noDatr(i), noDatc(i)) = vals(datr(nearestPtNdx(i)), datc(nearestPtNdx(i)));
        end
        
        % smooth the landscape
        pdN = 10;
        vals_pad = padarray(vals,[pdN pdN],'replicate','both');
        newState = conv2(vals_pad, ones(5,1)/5, 'same');
        newState = newState(pdN+1:end-pdN,pdN+1:end-pdN);
        
        newStateAll{transitions(n)} = newState;
        Vq{transitions(n)} = @(X,Y) interp2(XX,YY,newState,X,Y);
        %     [X,Y] = meshgrid([-5:0.1:5],[0:0.5:45]);
        %     figure;surf(X,Y,Vq(X,Y));view(45,45)%view(2);
        %     xlim([-5 5]);ylim([0 46]);zlim([0 1])
    end
else
    transitions =1:numel(newStateProbBaseline);
    fprintf('Not enough data to estimate Transition Probability.\n')
    fprintf('Setting Transition Probability to before first entry.\n')
    for n = 1:numel(newStateProbBaseline)
        newState = newStateProbBaseline(n).*ones(size(XX));
        newStateAll{transitions(n)} = newState;
        Vq{transitions(n)} = @(X,Y) interp2(XX,YY,newState,X,Y);
    end
end

end