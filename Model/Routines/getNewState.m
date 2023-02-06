function [newState] = getNewState(currState,TP,before,f,df,fBaseline)

nStateTransitions = numel(currState);
nStates = 3;
stateThresh = zeros(nStateTransitions,nStates);


% during
uniqueCurrStates = unique(currState(~before));
for state = 1:numel(uniqueCurrStates)
    currTPfunctions = TP.during(uniqueCurrStates(state),:);
    currFlyNdx = currState==uniqueCurrStates(state);
    nTransitions = numel(currTPfunctions);
    for newState = 1:nTransitions
        if ~isempty(currTPfunctions{newState})
            stateThresh(currFlyNdx,newState) = currTPfunctions{newState}(df(currFlyNdx),f(currFlyNdx));
        end
    end
end

% before
%nBaseline = sum(before);
stateThresh(before,:) = TP.before(currState(before),:);

% during baseline
if ~isempty(uniqueCurrStates)
    during_baseline = abs(df)<eps & abs(f-fBaseline)<eps;
    stateThresh(during_baseline,:) = TP.duringBaseline(currState(during_baseline),:);
end


stateThreshCum = cumsum(stateThresh,2)>=rand(nStateTransitions,1);
[~, newState] = max(stateThreshCum,[],2);

end