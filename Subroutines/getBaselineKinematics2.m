function [baseline] = getBaselineKinematics2(f_orco,border,condNdx,key)

% get the firing rate and change in firing rate
spk = f_orco.spk;
dSpk = f_orco.calcDeltaFR;

%close all
fe = f_orco.getFirstEntry('H',border);

if isempty(condNdx) && isempty(key)
    % separate out into each time point into before first entry (FE), below
    % baseline firing rate, baseline firing rate after FE, and above
    % baseline firing rate
    baselineSpk = f_orco.spk(1);
    condNdx = zeros(size(spk));condNdx(spk>baselineSpk) = 4;
    condNdx(spk<baselineSpk) = 2;condNdx(abs((spk-baselineSpk))<0.001 & dSpk==0) = 3;
    for j = 1:f_orco.nFly
        condNdx(j,1:fe(j)-1) = 1;
    end
    key = {'before','below','baseline','above'};
end

maxTimeDuring = f_orco.nPt-min(fe)+1;
maxTimeBefore = max(fe);

[condState, ~] = f_orco.getCond(condNdx,key);

tmp = cell(16,1);
tmp2 = cell(4,1);
for i = 1:f_orco.nFly
    [startNdx,endNdx,type] = startEndSeq(condState(i,:));
    for j = 1:max(type)
        tmp{j} = [tmp{j} endNdx(type == j)-startNdx(type == j)+1];
    end
    
    [startNdx,endNdx,type] = startEndSeq(f_orco.states.ndx(i,:));
    for j = 1:4
        tmp2{j} = [tmp2{j} endNdx(type == j)-startNdx(type == j)+1];
    end
end

fields = {'spk','curv','curvSmooth','spd','dur','state'};

for i = 1:numel(fields)
    baseline.before.([fields{i}]) = nan(f_orco.nFly,30,maxTimeBefore);
    baseline.during.([fields{i}]) = nan(f_orco.nFly,30,maxTimeDuring);
    baseline.first.([fields{i}]) = nan(f_orco.nFly,1,maxTimeDuring);
end

trackType = 'baseline';
% for before period
NdxType = find(cellfun(@(x) strcmpi(x,'before'), key));
[baseline.before] = getBaseline(f_orco,condNdx,NdxType,baseline.before,0);

% for last track of baseline in the during period
NdxType = find(cellfun(@(x) strcmpi(x,trackType), key));
[baseline.during] = getBaseline(f_orco,condNdx,NdxType,baseline.during,0);

% % for last track of baseline in the during period
% NdxType = find(cellfun(@(x) strcmpi(x,trackType), key));
% [baseline.first] = getBaselineAll(f_orco,condNdx,NdxType,baseline.first,fields,1);

end

function [baselineType] = getBaseline(f_orco,condNdx,NdxType,baselineType,n)
spd = f_orco.calcSpd;
phi = f_orco.calcPhi.*180./pi;
dur = f_orco.getCondDur./f_orco.fs;
[ang,dAng] = f_orco.calcAng;
%dAng = dAng.*180./pi;dAng = mod(dAng,360);
bordCond = f_orco.states.ndx==find(strcmpi(f_orco.states.key,'boundary'));
stopCond = f_orco.states.ndx==find(strcmpi(f_orco.states.key,'stops'));

for i = 1:f_orco.nFly
    [startNdx,endNdx,type] = startEndSeq(condNdx(i,:)==NdxType);
    startNdx = startNdx(type==1);endNdx = endNdx(type==1);
    len = endNdx-startNdx+1;
    startNdx(len<3) = [];endNdx(len<3) = [];
    len = endNdx-startNdx+1;
    
    if numel(startNdx)>n
        for j = 1:numel(startNdx)
            tmpCurv = (f_orco.curv(i,startNdx(j):endNdx(j))).*180./pi;%abs
            tmpPhi = (phi(i,startNdx(j):endNdx(j)));%abs
            tmpdAng = (dAng(i,startNdx(j):endNdx(j)));%abs

            tmpBordCond = bordCond(i,startNdx(j):endNdx(j));
            tmpStopCond = stopCond(i,startNdx(j):endNdx(j));
            tmpCurv(tmpBordCond) = tmpPhi(tmpBordCond);
            %tmpCurv(tmpStopCond) = tmpdAng(tmpStopCond);

            kernel = [ones(1,15)]./15;
            smoothCurv = conv(tmpCurv,kernel,'same');
            %figure;plot(tmpCurv);hold on;plot(smoothCurv);plot(smooth(tmpCurv,30))

            baselineType.spk(i,j,1:len(j)) = f_orco.spk(i,startNdx(j):endNdx(j));
            baselineType.curv(i,j,1:len(j)) = tmpCurv;
            baselineType.curvSmooth(i,j,1:len(j)) = smoothCurv;
            baselineType.spd(i,j,1:len(j)) = spd(i,startNdx(j):endNdx(j));
            baselineType.dur(i,j,1:len(j)) = dur(i,startNdx(j):endNdx(j));
            baselineType.state(i,j,1:len(j)) = f_orco.states.ndx(i,startNdx(j):endNdx(j));
        end
    end
end

end

function [baselineType] = getBaselineAll(f_orco,condNdx,NdxType,baselineType,fields,n)
spd = f_orco.calcSpd;
phi = f_orco.calcPhi.*180./pi;
dur = f_orco.getCondDur./f_orco.fs;
[ang,dAng] = f_orco.calcAng;
%dAng = dAng.*180./pi;dAng = mod(dAng,360);%(dAng>360) = dAng(dAng>360)-360;
bordCond = f_orco.states.ndx==find(strcmpi(f_orco.states.key,'boundary'));
stopCond = f_orco.states.ndx==find(strcmpi(f_orco.states.key,'stops'));

for i = 1:f_orco.nFly
    [startNdx,endNdx,type] = startEndSeq(condNdx(i,:)==NdxType);
    startNdx = startNdx(type==1);endNdx = endNdx(type==1);
    len = endNdx-startNdx+1;
    startNdx(len<3) = [];endNdx(len<3) = [];
    len = endNdx-startNdx+1;
    
    if numel(startNdx)>n
        len2 = endNdx(end)-startNdx(1)+1;
        allNdx = startNdx(1):endNdx(end);
        tmpCurv = (f_orco.curv(i,allNdx)).*180./pi;%abs
        tmpPhi = (phi(i,allNdx));%abs
        tmpdAng = (dAng(i,allNdx));%abs
        
        tmpBordCond = bordCond(i,allNdx);
        tmpStopCond = stopCond(i,allNdx);
        tmpCurv(tmpBordCond) = tmpPhi(tmpBordCond);
        %tmpCurv(tmpStopCond) = tmpdAng(tmpStopCond);% changed from dPhi to
        %change in sum of curvature
        
        
        kernel = [0 ones(1,15), 0]./15;
        smoothCurv = conv(tmpCurv,kernel,'same');
%         figure;plot(tmpCurv);hold on;plot(smoothCurv);plot(smooth(tmpCurv,30))

        baselineType.spk(i,1,1:len2) = f_orco.spk(i,allNdx);
        baselineType.curv(i,1,1:len2) = tmpCurv;
        baselineType.curvSmooth(i,1,1:len2) = smoothCurv;
        baselineType.spd(i,1,1:len2) = spd(i,allNdx);
        baselineType.dur(i,1,1:len2) = dur(i,allNdx);
        baselineType.state(i,1,1:len2) = f_orco.states.ndx(i,allNdx);
        
        GoodNdx = false(size(baselineType.spk(i,:)));
        for j = 1:numel(len)
            GoodNdx(startNdx(j)-startNdx(1)+1:endNdx(j)-startNdx(1)+1) = true;
        end
        
        for j = 1:numel(fields)
            try
            baselineType.([fields{j}])(i,1,~GoodNdx) = nan;
            catch
            end
        end
        
    end
end

end