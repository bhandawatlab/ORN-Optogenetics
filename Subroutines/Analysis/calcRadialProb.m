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
Before = cell(1,self.nFly);During = cell(1,self.nFly);
for i = 1:self.nFly
    BeforeTmp = self.(['r' lab])(i,1:fe(i));
    DuringTmp = self.(['r' lab])(i,fe(i):end);
    BeforeTmp(longStops(i,1:fe(i))) = [];
    DuringTmp(longStops(i,fe(i):end)) = [];
    
    Before{i} = BeforeTmp./self.rBound;
    During{i} = DuringTmp./self.rBound;
end
nBefore = cellfun(@(x) numel(x),Before)';
nDuring = cellfun(@(x) numel(x),During)';
w_Before = nBefore./sum(nBefore);
w_During = nBefore./sum(nDuring);

edges = 0:0.05:1;

% weighted standard dev and average across flies
N_byFly_bef = zeros(self.nFly,numel(edges)-1);
N_byFly_dur = zeros(self.nFly,numel(edges)-1);
for i = 1:self.nFly
    [N_byFly_bef(i,:),~] = histcounts(Before{i},edges);
    [N_byFly_dur(i,:),~] = histcounts(During{i},edges);
end
N_byFly_bef = N_byFly_bef+eps;N_byFly_bef = N_byFly_bef./sum(N_byFly_bef,2);
N_byFly_dur = N_byFly_dur+eps;N_byFly_dur = N_byFly_dur./sum(N_byFly_dur,2);
radialProb.before.flyMu = nanmean(N_byFly_bef);
radialProb.before.std = std(N_byFly_bef,w_Before);
%radialProb.before.std = nanstd(N_byFly_bef);
radialProb.before.sem = radialProb.before.std./(self.nFly-1);
radialProb.during.std = std(N_byFly_dur,w_During);
%radialProb.during.std = nanstd(N_byFly_dur);
radialProb.during.flyMu = nanmean(N_byFly_dur);
radialProb.during.sem = radialProb.during.std./(self.nFly-1);

% average across all data points (each fly is weighted by amount of time)
[N,~] = histcounts(cell2mat(Before),edges);
N = N+eps;
radialProb.before.yRaw = Before;
radialProb.before.weightedMu = N./sum(N);
radialProb.before.x = edges(2:end)-(edges(2)-edges(1))./2;
[N,~] = histcounts(cell2mat(During),edges);
N = N+eps;
radialProb.during.yRaw = During;
radialProb.during.weightedMu = N./sum(N);
radialProb.during.x = edges(2:end)-(edges(2)-edges(1))./2;

end