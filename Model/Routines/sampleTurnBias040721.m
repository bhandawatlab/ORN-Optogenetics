function [direction,opt] = sampleTurnBias040721(self,state,t,f,df,tt,xPos,yPos,...
    tSinceInhibition,sampCurv,sampSpeed,sampDur,sce)
%t = t.*0+1;
nFly = numel(t);
nPt = size(xPos,2);
m = self.model.TurnBias;

turnBias = zeros(size(t));
[~,tSlice]=max(t<(m.tt+m.tt(2)./2),[],2);
uSlice =unique(tSlice);
for s = 1:numel(uSlice)
    slice = uSlice(s);
    turnBias(tSlice==slice) = m.during{slice,state}(df(tSlice==slice),f(tSlice==slice));
end
turnBias(t<0) = m.before{state};
turnBias(tt>0) = m.during_baseline{state};%m.before{state};

inhibitionTBNdx = ~isnan(tSinceInhibition);
turnBias_inh = turnBias;
if any(inhibitionTBNdx) && ~isempty(m.inhibitionKin{state}{1})
    for s = 1:numel(uSlice)
        slice = uSlice(s);
        turnBias_inh(tSlice==slice) = m.inhibitionKin{state}{slice}(tSinceInhibition(tSlice==slice));
    end
%     inhibitionTB = m.inhibitionKin{state}(tSinceInhibition);%+1 = force optimal
%     turnBias(inhibitionTBNdx) = inhibitionTB(inhibitionTBNdx);
end
turnBias(inhibitionTBNdx) = turnBias_inh(inhibitionTBNdx);

[~,thetaAll] = calcCurv(xPos,yPos,nFly,nPt);
thetaAll = thetaAll+pi;
currAng = thetaAll;
currAng = nanmean(currAng,2);
currAng = wrapTo360(currAng.*180./pi)'-90;
%--------------------------------------------------------------------------
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

if any(isnan(turnBias))
    disp('meow walk')
    turnBias(isnan(turnBias)) = 1;
end
%------------------
%turnBias(isnan(turnBias)) = 1;
%------------------

% sample whether or not it's optimal
direction = opt.*sign((rand(nFly,1)<turnBias)-0.5);%- added negative?

%--------------------------------------------------------------------------
% % using turn bias rather than turn dir
% [direction] = getTurnBias([xPos(:,end),yPos(:,end)],currAng',sampCurv,turnBias,sampSpeed,sampDur,sce);
% direction = sign(direction);
%--------------------------------------------------------------------------



end

% [~, thetaAll] = self.calcCurv;
% l = length(s);
% currAng = thetaAll(fly,max(l-hist,1):l);
% currAng = nanmean(currAng);
% currAng = wrapTo360(currAng.*180./pi);
% prevDir = [cosd(currAng);sind(currAng)];
% dirRelCent = -[self.xH(fly,prevEndNdx); self.yH(fly,prevEndNdx)];
% n = [0;0;1];
% v = [prevDir;0];
% u = [dirRelCent;0];
% theta = get2VecAngle(v,u,n,false);

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





% %--------------------------------------------------------------------------
% x = xyPos(:,1); y = xyPos(:,2);
% nFlys = size(x,1);
% %u = [cosd(currAng),sind(currAng),zeros(nFlys,1)];
% v1 = [cosd(currAng+sampYaw),sind(currAng+sampYaw),zeros(nFlys,1)];
% v2 = [cosd(currAng-sampYaw),sind(currAng-sampYaw),zeros(nFlys,1)];
% w = [-x,-y,zeros(nFlys,1)];
% 
% % calculate the angles between the center facing vector and the
% % direction the fly is moving in before the turn and after the turn
% a = cross(v1,w);
% b = dot(v1',w')';
% c = sqrt(sum(a.^2,2));
% turnLeft = atan2(c,b).*180./pi;
% a = cross(v2,w);
% b = dot(v2',w')';
% c = sqrt(sum(a.^2,2));
% turnRight = atan2(c,b).*180./pi;
% 
% % if any(turnRight<0) || any(turnLeft<0)
% %     a
% % end
% 
% turnDir = sampYaw;
% turnDir(turnRight<turnLeft) = -sampYaw(turnRight<turnLeft);
% 
% [sigTemp] = computeDir(radPos,pTurnIn,turningBias);
% turnDir(sigTemp) = turnDir(sigTemp).*-1;
% %--------------------------------------------------------------------------

            
% if any(abs(f-4.71)>0.1)
%     for i = 1:numel(f)
%         figure(100);
%         plot(xPos(i,:),yPos(i,:),'b','Linewidth',2);hold on;
%         q0 = quiver(xPos(i,end),yPos(i,end),-xPos(i,end),-yPos(i,end));
%         q0.Color = 'blue';q0.MaxHeadSize = 1;q0.LineStyle = '--';
%         q1 = quiver(xPos(i,end),yPos(i,end),cosd(currAng(i))./2,sind(currAng(i))./2);
%         q1.Color = 'blue';q1.MaxHeadSize = 1;q1.LineStyle = '--';
%         plotCircle([0,0],1,100,'r');
%         plotCircle([0,0],0.01,100,'r');
%         axis([-2 2 -2 2]);axis square
%         hold off;
%         pause(0.1);
%     end
% end

% %--------------------------------------------------------------------------
% [~,thetaAll] = calcCurv(xPos,yPos,nFly,nPt);
% 
% if isempty(thetaAll)
%     thetaAll = zeros(nFly,1);
% end
% currAng = thetaAll;
% currAng = nanmean(currAng,2);
% currAng = wrapTo360(currAng.*180./pi)';
% prevDir = [cosd(currAng);sind(currAng)];
% dirRelCent = -[xPos(:,end),yPos(:,end)]';
% n = [0;0;1];
% v = [prevDir;zeros(1,nFly)];
% u = [dirRelCent;zeros(1,nFly)];
% opt = zeros(nFly,1);
% for fly = 1:nFly
%     theta = get2VecAngle(v(:,fly),u(:,fly),n,false);
%     opt(fly) = sign(nanmean(theta));
% end
% opt(opt==0) = 1;
% %--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% n = [0;0;1];
% opt = zeros(nFly,1);
% for fly = 1:nFly
%     prevDir = [xPos(fly,end)-xPos(fly,1);
%                 yPos(fly,end)-yPos(fly,1)];
%     dirRelCent = -[xPos(fly,end); yPos(fly,end)];
%     
%     v = [prevDir;0];
%     u = [dirRelCent;0];
%     % get the optimal direction
%     opt(fly) = sign(get2VecAngle(v,u,n,false));
% end
% opt(opt==0) = 1;
%--------------------------------------------------------------------------