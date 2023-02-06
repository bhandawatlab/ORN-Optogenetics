function [states, key] = getStates(curvPks,curvWalks,stopCond,boundCond,npt)
% getStates  takes in structures containing indexes for each movement state
%   and reformats it into a matrix with each state given a different index
%
%   Inputs: curvPks = structure with field all that contain the indexes of
%               when each sharp turn trajectory occurs
%           curvWalks = structure with field all that contain the indexes 
%               of when each curved walk trajectory occurs
%           stopCond = structure with field all that contain the indexes of
%               when each stop trajectory occurs
%           boundCond = structure with field all that contain the indexes 
%               of when each boundary trajectory occurs
%           npt = number of time points for each experiment
%
%   Output: states = fly x time matrix of numbers corresponding to which
%               state occurs at each time point for each fly
%           key = cell array stating what state each number in states
%               correspond to 
%

nFly = numel(curvPks.all);
states = zeros(nFly,npt);
for i = 1:nFly
    for j = 1:numel(curvPks.all{i})
        states(i,curvPks.all{i}{j}(:,2)) = 1;
    end
    for j = 1:numel(curvWalks.all{i})
        states(i,curvWalks.all{i}{j}(:,2)) = 2;
    end
    for j = 1:numel(stopCond.all{i})
        states(i,stopCond.all{i}{j}(:,2)) = 3;
    end
    for j = 1:numel(boundCond.all{i})
        states(i,boundCond.all{i}{j}(:,2)) = 4;
    end
    [startNdx,endNdx,type] = startEndSeq(states(i,:)>0);
    startNdx(type) = [];endNdx(type) = [];
    len = endNdx-startNdx+1;
    %------------------------------------------------------------------
    if ~isempty(startNdx)
        for j = 1:numel(startNdx)
            states(i,startNdx(j):floor(len./2)+startNdx(j)) = states(i,max(startNdx(j)-1,1));
            try
                states(i,floor(len./2)+startNdx(j):endNdx(j)) = states(i,endNdx(j)+1);
            catch
                states(i,floor(len./2)+startNdx(j):endNdx(j)) = states(i,startNdx(j)-1);
            end
        end
    end
    %------------------------------------------------------------------
end
states = states(:,1:npt);
states = [states,states(:,end)];
key = {'sharp turns';'curved walks';'stops';'boundary'};
end