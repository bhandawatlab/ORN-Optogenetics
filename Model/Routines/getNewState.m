function [newState] = getNewState(currState,TP,before,f,df)
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
nBaseline = sum(before);

stateThresh(before,:) = TP.before(currState(before),:);


stateThreshCum = cumsum(stateThresh,2)>=rand(nStateTransitions,1);
[~, newState] = max(stateThreshCum,[],2);

end