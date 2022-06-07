% Demo code for multiple subplots of estimation plots. Original toolbox
% written by Adam Claridge-Chang's Lab
% 
% Liangyu Tao June 10, 2019. Can now do multiple subplots

clear;close all
addpath(genpath(pwd))

% load data
load fisheriris
% create an extra category (to be used for the pair wise multiple
% comparison plots)
species2 = species;
species2(40:50) = {'tmp'};

% define basic paramters
lims = [];                      % y limits to set (example: [0 1])
isPaired = 'N';                 % equivalent to paired t-test
circleSize = 60;                % size of scatter plot points
barstate = 'off';                % 'On' = bar graph, 'Off' = scatter plot

% define subplots (#row, #column,subplot #, figure #)
subplots = [4,3,1,1];
xSubplots = repmat(100./subplots(1),1,subplots(1));
ySubplots = repmat(100./subplots(2),1,subplots(2));

% for comparisons between two items
%--------------------------------------------------------------------------
figure(1);set(gcf,'Position',[842 42 838 924])
p = panel();
p.pack(xSubplots, ySubplots);
pwmd = [];
for i = 1:12
    subplots(3) = i;
    [ss,p,pwmd] = dabest2(meas(1:100,1),species(1:100),p,pwmd,...
        lims,isPaired,circleSize,barstate,subplots);
    % write the 95% confidence bounds on the plot
    confBound = [num2str(ss.mdCi(1)) '; ' num2str(ss.md) '; ' num2str(ss.mdCi(2))];
    title(confBound)
    % set axis labels for comparisons against category 1
    [I,J] = ind2sub(subplots(1:2),subplots(3));
    p(I,J).select();
    ylabel('Sepal length');
    xtickangle(10)
    title('Iris Flower comparisons','Interpreter','none')
end
%--------------------------------------------------------------------------
%
% for comparisons between odd number of items
%--------------------------------------------------------------------------
figure(2);set(gcf,'Position',[842 42 838 924])
subplots(4) = 2;
p = panel();
p.pack(xSubplots, ySubplots);
pwmd = [];
for i = 1:12
    subplots(3) = i;
    [ss,p,pwmd] = dabest2(meas(:,1),species,p,pwmd,...
        lims,isPaired,circleSize,barstate,subplots);
    % set axis labels for comparisons against category 1
    [I,J] = ind2sub(subplots(1:2),subplots(3));
    p(I,J,1,1).select();
    ylabel('Sepal length');
    title('Iris Flower comparisons','Interpreter','none')
    p(I,J,2,1).select();
    ylabel('delta Sepal length');
    xtickangle(10)
end
suptitle('fisheriris')
%--------------------------------------------------------------------------
%
% for comparisons between even number of items >2
%--------------------------------------------------------------------------
figure(3);set(gcf,'Position',[842 42 838 924])
subplots(4) = 3;
p = panel();
p.pack(xSubplots, ySubplots);
% for pair wise multiple comparison
figure(4);set(gcf,'Position',[842 42 838 924])
pwmd = panel();
pwmd.pack(xSubplots, ySubplots);
for i = 1:12
    subplots(3) = i;
    [ss,p,pwmd] = dabest2(meas(:,1),species2,p,pwmd,...
        lims,isPaired,circleSize,barstate,subplots);
    % set axis labels for comparisons against category 1
    [I,J] = ind2sub(subplots(1:2),subplots(3));
    p(I,J,1,1).select();
    ylabel('Sepal length');
    title('Iris Flower comparisons','Interpreter','none')
    p(I,J,2,1).select();
    ylabel('delta Sepal length');
    xtickangle(10)
    % set axis labels for pair wise multiple comparison
    pwmd(I,J,1,1).select();
    ylabel('Sepal length');
    title('Iris Flower PWMC','Interpreter','none')
    pwmd(I,J,2,1).select();
    ylabel('delta Sepal length');
    xtickangle(10)
end
figure(3);
suptitle('fisheriris')
figure(4);
suptitle('fisheriris')
%--------------------------------------------------------------------------

