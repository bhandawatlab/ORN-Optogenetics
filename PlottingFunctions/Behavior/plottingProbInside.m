function [sliding_mean] = plottingProbInside(empFlys,border,numfliesEmp,fs,startPts,type)
if strcmpi(type,'during')
    [bin_inside_rim_sliding_SEM,bin_inside_rim_sliding_mean,inside_rim,framesinBin] = analysisFunc(empFlys,border,numfliesEmp,startPts,type);
    endPts = size(empFlys.rHead,2)-startPts+1;
    plottingFunc2(bin_inside_rim_sliding_SEM,bin_inside_rim_sliding_mean,inside_rim,framesinBin,fs,endPts)
else
    startPts = ones(size(startPts));
    [bin_inside_rim_sliding_SEM,bin_inside_rim_sliding_mean,inside_rim,framesinBin] = analysisFunc(empFlys,border,numfliesEmp,startPts,type);
    plottingFunc(bin_inside_rim_sliding_SEM,bin_inside_rim_sliding_mean,inside_rim,framesinBin,fs)
end

tmp = bin_inside_rim_sliding_mean';
 
len = min(length(tmp));
sliding_mean = tmp(1:len);
 
end
 
function [bin_inside_rim_sliding_SEM,bin_inside_rim_sliding_mean,inside_rim,framesinBin] = analysisFunc(Flys,border,numflies,startPts,type)
 
inside_rimTmp = Flys.rHead'<border;
time_to_end = size(inside_rimTmp,1)-min(startPts)+1;
inside_rim = nan(time_to_end,numflies);
for fly = 1:numflies
    tmp = inside_rimTmp(startPts(fly):end,fly);
    inside_rim(1:numel(tmp),fly) = tmp;
end
inside_rim_mean = nanmean(inside_rim,2);
framesinBin = 6; %6 frames = 200 ms
 
%==========================================================================
%sliding bins with 6 frames/bin
bin_inside_rim_sliding = zeros(length(inside_rim_mean),numflies);
 
for fly=1:numflies
    for i = 1:length(inside_rim_mean)
        if i< 3 % get the average of frames 1~i
        bin_inside_rim_sliding(i,fly) = mean(inside_rim(1:i,fly));  
        elseif i > (length(inside_rim_mean) -3)
            bin_inside_rim_sliding(i,fly) = nanmean(inside_rim(i-2:end,fly));
        else
            bin_inside_rim_sliding(i,fly) = nanmean(inside_rim(i-2:i+3,fly));
        end
    end
end
         
%mean of binned data
bin_inside_rim_sliding_mean =  nanmean(bin_inside_rim_sliding,2);
 
%std of binned data
bin_inside_rim_sliding_std = nanstd(bin_inside_rim_sliding,0,2);
 
%SEM
bin_inside_rim_sliding_SEM = bin_inside_rim_sliding_std./(sqrt(numflies));
 
end
 
function [] = plottingFunc(bin_inside_rim_sliding_SEM,bin_inside_rim_sliding_mean,inside_rim,framesinBin,fs)
%==========================================================================

%convert x range from frame # to time (min)
xrange_min = linspace(1/fs,(size(inside_rim,1)/(60*fs)),size(inside_rim,1));
 
SEM_y_plot = [bin_inside_rim_sliding_mean'- bin_inside_rim_sliding_SEM';(2*bin_inside_rim_sliding_SEM')];
h = area(xrange_min,SEM_y_plot');
set(h(1),'visible','off');
set(h(2),'FaceColor','k','EdgeColor','none');
alpha(.2);
hold on
 
plot(xrange_min,bin_inside_rim_sliding_mean,'k','linewidth',1.5);
 
plot([3 3],[0 1],'r:');%odor on
plot([6 6],[0 1],':'); %odor off
 
%xlabel('time (min)');
%ylabel('probability of being inside');
 
ylim([0 1]);
set(gca,'box','off','tickdir','out','ytick',[0:.2:1]);
%title(['Probability of being inside (sliding bins), ' num2str(framesinBin/fs*1000) ' msec/bin ' ]);
 
%set(gcf,'position',[34 518 819 420])
 
end

function [] = plottingFunc2(bin_inside_rim_sliding_SEM,bin_inside_rim_sliding_mean,inside_rim,framesinBin,fs,endPts)
%==========================================================================

%convert x range from frame # to time (min)
xrange_min = linspace(1/fs,(size(inside_rim,1)/(60*fs)),size(inside_rim,1));
 
SEM_y_plot = [bin_inside_rim_sliding_mean'- bin_inside_rim_sliding_SEM';(2*bin_inside_rim_sliding_SEM')];
h = area(xrange_min,SEM_y_plot');
set(h(1),'visible','off');
set(h(2),'FaceColor','k','EdgeColor','none');
alpha(.2);
hold on
 
plot(xrange_min-1/fs,bin_inside_rim_sliding_mean,'k','linewidth',1.5);
endPts = sort(endPts);
for i = 1:numel(endPts)
    plot([endPts(i) endPts(i)]./(60*fs),[0 1],':'); %odor off
    text(endPts(i)./(60*fs), 0.95, num2str(numel(endPts)-i))
end
% plot([3 3],[0 1],'r:');%odor on
% plot([6 6],[0 1],':'); %odor off
 
%xlabel('time (min)');
%ylabel('probability of being inside');
 
ylim([0 1]);
set(gca,'box','off','tickdir','out','ytick',[0:.2:1]);
%title(['Probability of being inside (sliding bins), ' num2str(framesinBin/fs*1000) ' msec/bin ' ]);
 
%set(gcf,'position',[34 518 819 420])
 
end