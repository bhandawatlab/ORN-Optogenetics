function [] = violinPlotsStats(data,opts)
% violinPlotsStats  significant box of violins :)
%
% Inputs:
%    data: A matrix of data. Each row is an observation. each column is a
%       variable.
%    opts: a structure of options with potential fields:
%       yl: y-limit (e.g. [ymin ymax])
%       sigLevel: significance level for sign rank or rank sum test
%       xLabels: cell array of labels for each variable
%       omitOutliers: true/(false) omitt outliers
%       isPaired: true/(false) whether each variable as paired
%       plotBox: (true)/false plot box plots
%       plotCentralTendency: true/(false) plot the median as horizontal bars
%       plotViolin: (true)/false plot the violin plot
%       plotData: (true)/false plot a scatter plot of the raw data points
%       plotStaggered: true/(false) stagger the x position of the raw data
%           points based on density estimation
%       detault for boolean options is in ()
%
% Example usage:
%       x = rand(100,2);
%       violinPlotsStats(x,[]);
%

% setup commonly used values
nCond = size(data,2);
nObs = size(data,1);
nPt = 200;

if ~isfield(opts,'yl') || isempty(opts.yl)
    yl = [min(data(:)) max(data(:))];
    opts.yl = yl+(diff(yl)+1).*[-0.1 0.1];
end
if ~isfield(opts,'sigLevel') || isempty(opts.sigLevel)
    opts.sigLevel = [0.05,0.01,0.001];
end
if ~isfield(opts,'xLabels') || isempty(opts.xLabels)
    opts.xLabels = cellstr(num2str([1:nCond]'));
end
if ~isfield(opts,'omitOutliers') || isempty(opts.omitOutliers)
    opts.omitOutliers = false;
end
if ~isfield(opts,'plotBox') || isempty(opts.plotBox)
    opts.plotBox = true;
end
if ~isfield(opts,'plotCentralTendency') || isempty(opts.plotCentralTendency)
    opts.plotCentralTendency = false;
end
if ~isfield(opts,'plotViolin') || isempty(opts.plotViolin)
    opts.plotViolin = true;
end
if ~isfield(opts,'plotData') || isempty(opts.plotData)
    opts.plotData = true;
end
if ~isfield(opts,'isPaired') || isempty(opts.plotData)
    opts.plotPaired = false;
    opts.isPaired = false;
end
if ~isfield(opts,'plotStaggered') || isempty(opts.plotData)
    opts.plotStaggered = false;
end

cc = distinguishable_colors(nCond);

% statistical test
p = nan(nCond,nCond);
p2 = nan(nCond,nCond);
for i = 1:nCond
    for j = i+1:nCond
        p(i,j) = signrank(data(:,i),data(:,j));
        p2(i,j) = ranksum(data(:,i),data(:,j));
    end
end
if any(sum(isnan(data))./size(data,1)>0.5)
    p = p2;disp('Note enough paired data. Using rank sum')
    test = 'Rank Sum';
elseif ~opts.isPaired
    p = p2;
    test = 'Rank Sum';
else
    test = 'Sign Rank';
end
p(p>max(opts.sigLevel)) = nan;

if opts.omitOutliers
    outliers = isoutlier(data);
else
    outliers = false(size(data));
end


% get the density and randomized x point positions
% initialize matrix
f_dens = zeros(nCond,nPt);
x_dens = zeros(nCond,nPt);
x = zeros(nObs,nCond);
for c = 1:nCond
    currData = data(:,c);
    currData(outliers(:,c)) = nan;
    [f,~] = ksdensity(currData,data(:,c));
    [f_dens(c,:),x_dens(c,:)] = ksdensity(currData,linspace(opts.yl(1),opts.yl(2),200));
    %[f,~] = ksdensity(data(:,c),data(:,c),'Support','positive');
    %[f_dens,x_dens] = ksdensity(data(:,c),linspace(yl(1),yl(2),200),'Support','positive');
    m_f = max(f_dens(c,:));
    f = f./m_f.*0.8;
    if opts.plotStaggered
        x(:,c) = (rand(nObs,1)-0.5).*f+c;
    else
        x(:,c) = c;
    end
end

%% plotting
hold on;

% plot the data
if opts.plotData
    try
        scatter(x,data,'Color',cc,'LineWidth',1);%,'k'
    catch
        scatter(x,data,[],cc,'LineWidth',1);%,'k'
    end
end
if opts.plotPaired
    plot(x',data','k-','LineWidth',1);%,'k'
    
end

% plot the summary lines
if opts.plotCentralTendency
    %plot((1:4)+[-0.5; 0.5],repmat(nanmean(data),2,1),'k');
    plot((1:4)+[-0.5; 0.5],repmat(nanmedian(data),2,1),'--k');
end

if opts.plotBox
    boxplot(data,'Color',cc)
    %boxplot([1:nCond]+zeros(nObs,1),data)
end

% plot the violin plots
if opts.plotViolin
    m_f = max(f_dens,[],2);
    for c = 1:nCond
        tmp_f_dens = f_dens(c,:);
        tmp_f_dens(tmp_f_dens<0.001) = nan;
        plot([tmp_f_dens;-tmp_f_dens]./m_f(c)./2+c,x_dens(c,:),'-b','linewidth',1)
    end
end


% significance bar takes up 20% of the space
sigBar_lb = max(x_dens(f_dens>0.001));
sigBar_ub = (sigBar_lb-opts.yl(1)).*1.2;

[r,c] = find(~isnan(p));
d = c-r;
[~,ndx] = sort(d,'descend');
r = r(ndx);
c = c(ndx);
nSig = numel(c);
barOffset = (sigBar_ub-sigBar_lb)./nSig;
for i = 1:nSig
    x = [r(i) c(i)];
    y = sigBar_lb+barOffset.*(nSig-i);
    nDot = find(p(r(i),c(i))<=opts.sigLevel,1,"last");
    x_star = mean(x)+(0:nDot-1).*0.1-mean((0:nDot-1).*0.1);
    y_star = repmat(y+barOffset./3,1,nDot);
    plot(x,y.*[1 1],'k');
    plot(x_star,y_star,'k*')
    title(test)
end

xlim([0.5 nCond+0.5]);
ylim([opts.yl(1) max([sigBar_ub,opts.yl(2)])])
xticklabels(opts.xLabels)
hold off
end





