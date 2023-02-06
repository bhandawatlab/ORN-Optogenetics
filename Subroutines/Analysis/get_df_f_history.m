function [dFSmooth,FSmooth,ndx] = get_df_f_history(f_orco,dSpkRateAll,timeInterval)
% get_df_f_history  calculates the average firing rate over a period of
%   time.
%
%   Inputs: self = fly object
%           dSpkRateAll = change in firing rate
%           timeInterval = time inverval in ms to average over
%
%   Output: dFSmooth = average firing rate
%           FSmooth = average change in firing rate
%           ndx = ~transposed hankel matrix (index for smoothing)
%   

hist_start = ceil(timeInterval(1)./1000.*f_orco.fs);%max(1,
hist_end = ceil(timeInterval(2)./1000.*f_orco.fs)-1;
%dt = hist_end-hist_start+1;

ndx = [1:size(dSpkRateAll,2)]-[hist_start:hist_end]'-1;
ndx(:,sum(ndx<1)>0) = [];

%ndx = (hist_start:size(dSpkRateAll,2)-dt)+[0:dt]'-hist_end;
for fly = 1:f_orco.nFly
    curr_df = dSpkRateAll(fly,:);
    curr_f = f_orco.spk(fly,:);
    
    dFSmooth(fly,:) = mean(curr_df(ndx),1);
    FSmooth(fly,:) = mean(curr_f(ndx),1);
end
dFSmooth = [zeros(size(dFSmooth,1),hist_end+1),dFSmooth];
FSmooth = [zeros(size(FSmooth,1),hist_end+1),FSmooth];

% %dFSmooth = [zeros(size(dFSmooth,1),hist_end),dFSmooth]./(hist_end-hist_start+1);% average
% %FSmooth = [zeros(size(FSmooth,1),hist_end),FSmooth]./(hist_end-hist_start+1);% average

end
