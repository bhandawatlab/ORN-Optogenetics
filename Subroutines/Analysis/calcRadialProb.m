function radialProb = calcRadialProb(self,lab,border,stopSpd,dur)
% calcRadialProb  calculates the radial occupancy.
%
%   Inputs: self = fly object
%           lab = label ('' = body, 'H' = head)
%           border = light arena border (in cm)
%           stopSpd = speed threshold to consider as stops (in mm/s)
%           dur = duration that is considered as long stops (in sec)
%
%   Output: radialProb = structure for radial occupancy. x are the center
%       of each bin, y is the probability, and yRaw is the raw positions
%   

fe = getFirstEntry(self,lab,border);
stops = self.calcSpd<stopSpd;
longStops = false(size(stops));
for i = 1:self.nFly
    [startNdx,endNdx,type] = startEndSeq(stops(i,:));
    len = endNdx-startNdx+1;
    startNdx(len<dur | type == false) = [];
    endNdx(len<dur | type == false) = [];
    for j = 1:numel(startNdx)
        longStops(i,startNdx(j):endNdx(j)) = true;
    end
    %--------
    [startNdx,endNdx,type] = startEndSeq(longStops(i,:));
    startNdx(type==0) = [];endNdx(type==0) = [];
    assert(all(endNdx-startNdx+1>=dur));
    %--------
end

% radial probability
Before = [];During = [];
for i = 1:self.nFly
    BeforeTmp = self.(['r' lab])(i,1:fe(i));
    DuringTmp = self.(['r' lab])(i,fe(i):end);
    BeforeTmp(longStops(i,1:fe(i))) = [];
    DuringTmp(longStops(i,fe(i):end)) = [];
    
    Before = [Before BeforeTmp];
    During = [During DuringTmp];
end
Before = Before./4;During = During./4;
edges = 0:0.05:1;
[N,~] = histcounts(Before,edges);
N = N+eps;
radialProb.before.yRaw = Before;
radialProb.before.y = N./sum(N);
radialProb.before.x = edges(2:end)-0.025;
[N,~] = histcounts(During,edges);
N = N+eps;
radialProb.during.yRaw = During;
radialProb.during.y = N./sum(N);
radialProb.during.x = edges(2:end)-0.025;
end