function [direction] = sampleTurnBiasStops(self,state,t,f,df,tt,xPos,yPos,curv)
%t = t.*0+1;
nFly = numel(t);
nPt = size(xPos,2);
m = self.model.TurnBias;
% xGrid = [-5:0.01:5];yGrid = [0:0.1:50];[XX,YY] = meshgrid(xGrid,yGrid);
% mu = mean(m.during{state}(XX,YY),'all');
%turnBias = (m.during{state}(df,f)-0.5).*2+0.5;
turnBias = m.during{state}(df,f);%+0.1;%+0.2;%0.7.*ones(size(df));
%turnBias = 1-m.during{state}(df,f);
turnBias(t<0) = m.before{state};

turnBias(tt>0) = m.during_baseline{state};%m.before{state};

%--------------------------------------------------------------------------
% startAng = myatan((xPos(:,end)-xPos(:,end-1))',(yPos(:,end)-yPos(:,end-1))','degrees',2);
% endAng = startAng+curv';
% startDir = [cosd(startAng);sind(startAng)];
% endDir = [cosd(endAng);sind(endAng)];
% startDirRelCent = -[xPos(:,end), yPos(:,end)]';
% endDirRelCent = -[xPos(:,end), yPos(:,end)]';
% 
% n = [0;0;1];
% n = repmat(n,1,nFly);
% v = [startDir;zeros(1,nFly)];
% u = [startDirRelCent;zeros(1,nFly)];
% thetaStart = get2VecAngle(v,u,n,false);
% v = [endDir;zeros(1,nFly)];
% u = [endDirRelCent;zeros(1,nFly)];
% thetaEnd = get2VecAngle(v,u,n,false);
% thetaDiff = thetaEnd-thetaStart;% negative means turned more inwards
% opt = sign(thetaDiff)';
% opt(opt==0) = 1;

[~,thetaAll] = calcCurv(xPos,yPos,nFly,nPt);
thetaAll = thetaAll+pi;
currAng = thetaAll;
currAng = nanmean(currAng,2);
currAng = wrapTo360(currAng.*180./pi)'-90;
prevDir = [cosd(currAng);sind(currAng)];
dirRelCent = -[xPos(:,end),yPos(:,end)]';
n = [0;0;1];
v = [prevDir;zeros(1,nFly)];
u = [dirRelCent;zeros(1,nFly)];
opt = zeros(nFly,1);
for fly = 1:nFly
    theta = get2VecAngle(v(:,fly),u(:,fly),n,false);
    opt(fly) = sign(nanmean(theta));
end
opt(opt==0) = 1;
%--------------------------------------------------------------------------

% sample whether or not it's optimal
direction = opt.*sign((rand(nFly,1)<turnBias)-0.5);% added negative?


end

function [curv,thetaAll] = calcCurv(x,y,nFly,nPt)
curv = zeros(nFly,nPt-1);
thetaAll = zeros(nFly,nPt-1);
for i = 1:nFly
    % calculating curvature
    Vertices = horzcat(x(i,:)',y(i,:)');
    %get normal vectors
    N=LineNormals2D(Vertices);
    
    theta = zeros(1,nPt-1); %make it zero vector
    for p=1:(nPt-1)
        %because normal is a unit vector. The x-component can be used
        % to determine the angle it makes with x-axis.
        theta(p)=acos(N(p,1));
        % the if loop converts from 0 to pi to -pi to pi
        if N(p,1)>0 && N(p,2)<0
            theta(p)=-theta(p);
        elseif N(p,1)<0 && N(p,2)<0
            theta(p)=-theta(p);
        end
    end
    theta(isnan(theta)) = 0;
    theta = unwrap(theta);
    thetaAll(i,1:end) = theta;
    curv(i,1:end-1) = diff(theta);
end
curv = [curv, zeros(nFly,1)];
end