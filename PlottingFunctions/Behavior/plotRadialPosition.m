function [xCent,yCent,distByTime,perAtBound,radialOccupancy,r,tt] = ...
    plotRadialPosition(self,border,cond,plotFig,varargin)
% default values
dt = 2.*self.fs;
dr = 0.2;
fe = self.getFirstEntry('H',border);
tMax = self.nPt./2;
clims = [0 1];
stopThresh = 0;
xrange = [0 tMax];yrange = [0 4];nOvlapX = 2*dt/4;nOvlapY = 2*dr/4;%2*dx/4;2*dy/4; 
plotType = 'surf';%'imagesc'

% user variable input
if numel(varargin)>0
    assert(mod(numel(varargin),2)==0)
    for i = 1:2:numel(varargin)-1
        if strcmpi(varargin{i},'dt')
            dt = varargin{i+1};
        elseif strcmpi(varargin{i},'dr')
            dr = varargin{i+1};
        elseif strcmpi(varargin{i},'fe')
            fe = varargin{i+1};
        elseif strcmpi(varargin{i},'timeBound')
            tMax = varargin{i+1};
        elseif strcmpi(varargin{i},'clims')
            clims = varargin{i+1};
        elseif strcmpi(varargin{i},'stopThresh')
            stopThresh = varargin{i+1};
        elseif strcmpi(varargin{i},'plotType')
            plotType = varargin{i+1};
        end
    end
end

% align by first entry
if strcmpi(cond,'before')
    r = nan(self.nFly,min(fe)-1);
    speed = self.calcSpd;
    for fly = 1:self.nFly
        tmp = self.rH(fly,1:fe(fly)-1);
        r(fly,:) = tmp(1:min(fe)-1);
        spd(fly,1:min(fe)-1) = speed(fly,1:min(fe)-1);
    end
    spd = spd(:,1:min(tMax,size(r,2)));
    r = r(:,1:min(tMax,size(r,2)));
    r = r./max(max(r,[],2),self.rBound).*self.rBound;
    r(spd<stopThresh) = nan;
    %r(r>self.rBound) = self.rBound;
elseif strcmpi(cond,'after')
    r = nan(self.nFly,self.nPt-min(fe)+1);
    speed = self.calcSpd;
    for fly = 1:self.nFly
        tmp = self.rH(fly,fe(fly):end);
        r(fly,1:numel(tmp)) = tmp;
        spd(fly,1:numel(tmp)) = speed(fly,fe(fly):end);
    end
    spd = spd(:,1:min(tMax,size(r,2)));
    r = r(:,1:min(tMax,size(r,2)));
    r = r./max(max(r,[],2),self.rBound).*self.rBound;
    r(spd<stopThresh) = nan;
    %r(r>self.rBound) = self.rBound;
end

[y,xCent,yCent,xRange,yRange,raw] = slidingBinsHistcounts(repmat(1:1:size(r,2),size(r,1),1),...
    r.^2,dt,dr.*4,xrange,yrange.^2,nOvlapX,nOvlapY*4,false);%true


distByTime = zeros(1,size(xRange,2)-1);
perAtBound = zeros(1,size(xRange,2)-1);
%figure;
for i = 1:size(xRange,2)
    tmp = r(:,xRange(1,i)+1:min(xRange(2,i),size(r,2)));
    tmpSpd = spd(:,xRange(1,i)+1:min(xRange(2,i),size(r,2)));
    perAtBound(i) = sum(tmp>3.85,'all')./sum(~isnan(tmp),'all');
    
    distByTime2(i) = nanmean(tmp,'all');
    tmp(tmp>3.85) = nan;
    %tmp(tmpSpd<0.5) = nan;
    distByTime(i) = nanmean(tmp,'all');
    %histogram(tmp,sqrt(0:0.2:16),'Normalization','probability');
end

radialOccupancy = zeros(1,size(yRange,2)-1);
for i = 1:size(yRange,2)
    radialOccupancy(i) = sum(r.^2<=yRange(2,i) & r.^2>yRange(1,i),'all');
end

% not radial density
[x1,x2] = meshgrid(xCent./self.fs, yCent);
xi = [x1(:) x2(:)./4];
tt = repmat(1:1:size(r,2),size(r,1),1)./self.fs;

X = [tt(:),(r(:))];
X(any(isnan(X),2),:) = [];

% downsample to 200,000 data points
if size(X,1)>200000
    X = X(randsample(size(X,1),200000),:);
end
if strcmpi(plotType,'imagesc')
    N = ksdensity(X,xi);
    N = reshape(N,size(x1));
    if plotFig
        imagesc(xCent./self.fs, yCent./4,N,clims);colorbar
    end
else
    if plotFig
        ksdensity(X,xi);colorbar
        colormap(jet)
        view(45,45);
        set(gca,'Xdir','normal','Ydir','reverse');
        zlabel('Density')
        xlim([0 180]);ylim([0 4]);zlim(clims);
        xticks([0:60:180]); xticklabels(cellstr(num2str([0:60:180]')))
        yticks([0:2:4]); yticklabels(cellstr(num2str([0:2:4]')))
    end
end
xlabel('time')
ylabel('r')
shading interp

end
