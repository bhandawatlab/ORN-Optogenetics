function [turnDir] = getTurnBias(xyPos,currAng,sampCurv,pTurnIn,sampSpeed,sampDur,sce)
x = xyPos(:,1); y = xyPos(:,2);
nFlys = size(x,1);
switch sce
    case 'ST'
        % ST
        %u = [cosd(currAng),sind(currAng),zeros(nFlys,1)];
        v1 = [cosd(currAng+sampCurv),sind(currAng+sampCurv),zeros(nFlys,1)];
        v2 = [cosd(currAng-sampCurv),sind(currAng-sampCurv),zeros(nFlys,1)];
        w = [-x,-y,zeros(nFlys,1)];
        
        % calculate the angles between the center facing vector and the
        % direction the fly is moving in before the turn and after the turn
        a = cross(v1,w);
        b = dot(v1',w')';
        c = sqrt(sum(a.^2,2));
        turnLeft = atan2(c,b).*180./pi;
        a = cross(v2,w);
        b = dot(v2',w')';
        c = sqrt(sum(a.^2,2));
        turnRight = atan2(c,b).*180./pi;
    case 'CW'
        % CW
        xEnd = zeros(nFlys,2);yEnd = zeros(nFlys,2);
        for j = 1:nFlys
            xEnd(j,1) = x(j)+sum(sampSpeed(j).*cosd(currAng(j)+sampCurv(j).*[1:1:sampDur(j)]));
            yEnd(j,1) = y(j)+sum(sampSpeed(j).*cosd(currAng(j)+sampCurv(j).*[1:1:sampDur(j)]));
            xEnd(j,2) = x(j)+sum(sampSpeed(j).*cosd(currAng(j)-sampCurv(j).*[1:1:sampDur(j)]));
            yEnd(j,2) = y(j)+sum(sampSpeed(j).*cosd(currAng(j)-sampCurv(j).*[1:1:sampDur(j)]));
        end
        
        turnLeft = sqrt((xEnd(:,1)).^2+(yEnd(:,1)).^2);
        turnRight = sqrt((xEnd(:,2)).^2+(yEnd(:,2)).^2);
end

turnDir = sampCurv;
turnDir(turnRight<turnLeft) = -sampCurv(turnRight<turnLeft);
pSamp = rand(length(turnDir),1);
sigTemp = pSamp>=pTurnIn;
turnDir(sigTemp) = turnDir(sigTemp).*-1;

end

function [curvDir] = walkTurnBias(xyPos,currAng,sampCurv,sampSpeed,sampDur,radPos,pTurnIn,Flys,turningBias)
x = xyPos(:,1); y = xyPos(:,2);
nFlys = length(Flys);
xEnd = zeros(nFlys,2);yEnd = zeros(nFlys,2);
for f = 1:nFlys
    j = Flys(f);
    xEnd(f,1) = x(j)+sum(sampSpeed(j).*cosd(currAng(j)+sampCurv(j).*[1:1:sampDur(j)]));
    yEnd(f,1) = y(j)+sum(sampSpeed(j).*cosd(currAng(j)+sampCurv(j).*[1:1:sampDur(j)]));
    xEnd(f,2) = x(j)+sum(sampSpeed(j).*cosd(currAng(j)-sampCurv(j).*[1:1:sampDur(j)]));
    yEnd(f,2) = y(j)+sum(sampSpeed(j).*cosd(currAng(j)-sampCurv(j).*[1:1:sampDur(j)]));
end

turnLeft = sqrt((xEnd(:,1)).^2+(yEnd(:,1)).^2);
turnRight = sqrt((xEnd(:,2)).^2+(yEnd(:,2)).^2);

curvDir = sampCurv(Flys);
curvDir(turnRight<turnLeft) = -sampCurv(turnRight<turnLeft);

[sigTemp] = computeDir(radPos,pTurnIn,turningBias);

curvDir(sigTemp) = curvDir(sigTemp).*-1;

end


function [dir] = computeDir(rPos,pTurnIn)

sBLim = (0:0.2:size(pTurnIn,2)*0.2)./4;sBLim(end) = sBLim(end)+0.1;
if sum((rPos<sBLim(end) & rPos>=sBLim(1)))>0
    pSamp = rand(length(rPos),1);
    dir = false(size(pSamp));
    if all(pTurnIn == pTurnIn(1,1))
        dir = pSamp>=pTurnIn(1,1);
    else
        for j = 1:size(pTurnIn,2)
            bordLoc = (rPos<sBLim(j+1) & rPos>=sBLim(j));
            SampChoice = pSamp>pTurnIn(1,j);
            SampChoice2 = pSamp>pTurnIn(2,j);
            SampChoice3 = pSamp>0.5;
            
            assert(length(bordLoc)==length(SampChoice),'err1')
            assert(length(bordLoc)==length(turningBias),'err2')
        end
    end
end
end