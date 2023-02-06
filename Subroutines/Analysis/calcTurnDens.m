function turnDens = calcTurnDens(self,lab,s,key)
% calcTurnDens  calculates the turn density
%
%   Inputs: self = fly object
%           lab = label ('' = body, 'H' = head)
%           s = structure with a cell array of time when turns occur.
%           key = cell array key for s indicating the index of each state
%           Note that both s and key are gotten from:
%           self.getKinematicsCond(condNdx,key,thresh,timeInterval,varargin)
%
%   Output: turnDens = structure for turn density. x are the center
%       of each bin, Bef/Dur is the density, and BefRaw/durRaw are the raw
%       radial positions of turns
%   

dx = [0:0.05:1];%0:0.025:1;%
dx(end) = 1.1;

STNdx = any(strcmpi('sharp turns', key),2);
turns = s.ndx(:,STNdx);
cond = key(STNdx,2);nCond = numel(cond);
befNdx = strcmpi(cond,'before');

tmpBef = [];tmpDur = [];
for fly = 1:self.nFly
    tmpBef = [tmpBef self.(['r' lab])(fly,turns{fly,befNdx})];
    tmpDur = [tmpDur self.(['r' lab])(fly,cell2mat(reshape(turns(fly,~befNdx),[],1)))];
end
tmpBef = tmpBef./self.rBound;tmpDur = tmpDur./self.rBound;
decisionDensBef = histcounts(tmpBef,dx)./(dx(2:end).^2 - dx(1:end-1).^2);
turnDens.Bef = decisionDensBef./sum(decisionDensBef);
decisionDensDur = histcounts(tmpDur,dx)./(dx(2:end).^2 - dx(1:end-1).^2);
turnDens.Dur = decisionDensDur./sum(decisionDensDur);
turnDens.x = dx;

turnDens.BefRaw = tmpBef;
turnDens.durRaw = tmpDur;

end