function [] = plotDistributions(self,stateDistBefore,stateDistEarly,stateDistLate,earlyNdx,lateNdx,nMin,fNum)
% define basic paramters
lims = [];                      % y limits to set (example: [0 1])
isPaired = 'N';                 % equivalent to paired t-test
circleSize = 60;                % size of scatter plot points
barstate = 'off';                % 'On' = bar graph, 'Off' = scatter plot
subplots = [4,3,1,fNum]; % ysbplt, xsbplt, sbpltN, fig
xSubplots = repmat(100./subplots(1),1,subplots(1));
ySubplots = repmat(100./subplots(2),1,subplots(2));

figure(subplots(4));set(gcf,'Position',[842 42 838 924])
p = panel();
p.pack(xSubplots, ySubplots);


% plotting ST
key = {'before',['first ' num2str(earlyNdx(end))],'later 3-10'};
for i = 1:numel(stateDistBefore)
    stCurv = cell(0);stDur = cell(0);stSpd = cell(0);
    
    if isempty(stateDistEarly{i}.curv)
        stateDistEarly{i}.curv = zeros(1,50)-10.*eps.*rand(1,50);
    end
    if isempty(stateDistEarly{i}.dur)
        stateDistEarly{i}.dur = zeros(1,50)-10.*eps.*rand(1,50);
    end
    if isempty(stateDistEarly{i}.spd)
        stateDistEarly{i}.spd = zeros(1,50)-10.*eps.*rand(1,50);
    end
    if isempty(stateDistLate{i}.curv)
        stateDistLate{i}.curv = zeros(1,50)-10.*eps.*rand(1,50);
    end
    if isempty(stateDistLate{i}.dur)
        stateDistLate{i}.dur = zeros(1,50)-10.*eps.*rand(1,50);
    end
    if isempty(stateDistLate{i}.spd)
        stateDistLate{i}.spd = zeros(1,50)-10.*eps.*rand(1,50);
    end
    
    stCurv(1:numel(stateDistBefore{i}.curv),1) = num2cell(stateDistBefore{i}.curv);
    stCurv(1:numel(stateDistEarly{i}.curv),2) = num2cell(stateDistEarly{i}.curv);
    stCurv(1:numel(stateDistLate{i}.curv),3) = num2cell(stateDistLate{i}.curv);
    stDur(1:numel(stateDistBefore{i}.dur),1) = num2cell(stateDistBefore{i}.dur);
    stDur(1:numel(stateDistEarly{i}.dur),2) = num2cell(stateDistEarly{i}.dur);
    stDur(1:numel(stateDistLate{i}.dur),3) = num2cell(stateDistLate{i}.dur);
    stSpd(1:numel(stateDistBefore{i}.spd),1) = num2cell(stateDistBefore{i}.spd);
    stSpd(1:numel(stateDistEarly{i}.spd),2) = num2cell(stateDistEarly{i}.spd);
    stSpd(1:numel(stateDistLate{i}.spd),3) = num2cell(stateDistLate{i}.spd);

%     [subplots,p] = plottingFunc(self,stCurv,stDur,stSpd,key,p,lims,isPaired,...
%         circleSize,barstate,subplots,self.states.key{i})
    [subplots,p] = plottingFunc(self,stCurv,stDur,stSpd,key,p,lims,isPaired,...
        circleSize,barstate,subplots,self.states.key{i},nMin);
end

end



function [subplots,p] = plottingFunc(self,stCurv,stDur,stSpd,keyCond,p,lims,isPaired,...
    circleSize,barstate,subplots,titl,nMin)

nCond = size(stCurv,2);
for i = 1:nCond
    if numel(cell2mat(stCurv(:,i)))>=3
        stCurvAll{i,1} = cell2mat(stCurv(:,i));
        stDurAll{i,1} = cell2mat(stDur(:,i));
        stSpdAll{i,1} = cell2mat(stSpd(:,i));
        gCurv{i,1} = repmat(keyCond(i),numel(stCurvAll{i,1}),1);
        gDur{i,1} = repmat(keyCond(i),numel(stDurAll{i,1}),1);
        gSpd{i,1} = repmat(keyCond(i),numel(stSpdAll{i,1}),1);
    end
end

dat = abs(cell2mat(stCurvAll));g = cat(1,gCurv{:});
dat = mod(dat,360);
if strcmpi(titl,'sharp turns')
    currTitle = 'tot curv (deg)';
elseif strcmpi(titl,'curved walks')
    currTitle = 'avg curv (deg/s)';
elseif strcmpi(titl,'boundary')
    currTitle = 'reorient (deg)';
elseif strcmpi(titl,'stops')
    currTitle = 'phi (deg)';
end
l = [0 360];

[ss,p,~] = dabest3(dat+eps*rand(size(dat)),g,p,[],l,isPaired,circleSize,barstate,subplots);
labelAxis(self,p,subplots,currTitle,titl)
subplots(3) = subplots(3)+1;

if strcmpi(titl,'sharp turns')
    currTitle = 'tot curv (deg)';
elseif strcmpi(titl,'curved walks')
    currTitle = 'avg curv (deg/s)';
elseif strcmpi(titl,'boundary')
    currTitle = 'reorient (deg)';
elseif strcmpi(titl,'stops')
    currTitle = 'phi (deg)';
end

dat = cell2mat(stDurAll);g = cat(1,gDur{:});
[ss,p,~] = dabest3(dat,g,p,[],lims,isPaired,circleSize,barstate,subplots);
labelAxis(self,p,subplots,'dur (s)',titl)
subplots(3) = subplots(3)+1;

dat = cell2mat(stSpdAll);g = cat(1,gSpd{:});
[ss,p,~] = dabest3(dat,g,p,[],[0 30],isPaired,circleSize,barstate,subplots);
labelAxis(self,p,subplots,'spd (mm/s)',titl)
subplots(3) = subplots(3)+1;

end


function [] = labelAxis(self,p,subplots,label,titl)
[J,I] = ind2sub(subplots(2:-1:1),subplots(3));
try
    p(I,J).select();
    ylabel(label);
    title(titl,'Interpreter','none')
catch
    p(I,J,1,1).select();
    ylabel(label);
    title(titl,'Interpreter','none')
    p(I,J,2,1).select();
    ylabel(['delta ' label]);
    xtickangle(10)
end

end





