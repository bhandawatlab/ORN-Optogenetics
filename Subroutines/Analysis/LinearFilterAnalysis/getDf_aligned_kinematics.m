function [df_leaving_all,df_entering_all,f_leaving_all,f_entering_all,...
    kinematics_leaving,kinematics_entering,kinematics_leaving_std,kinematics_entering_std] = ...
    getDf_aligned_kinematics(f_orco,thresh,delay,yl,type,singleTrial)
%getDf_aligned_kinematics  calculates the firing rate, change in firing
%   rate, and the kinematics when leaving or entering (aligned by peak in
%   change in firing rate)
%
%   Inputs: f_orco = fly object
%           thresh = threshold for peak (or trough) in delta firing rate
%           delay = time window to consider before and after the peak
%           yl = ylabel (for plotting purposes)
%           type = string indicating to calculate curvature or speed 
%               ('curv' or 'spd')
%           singleTrial = true/false of whether to consider individual
%               flies separately
%
%   Output: df_leaving_all = vector of the average change in firing rate
%               at leaving (-delay to +delay)
%           df_entering_all = vector of the average change in firing rate
%               at entering
%           f_leaving_all = vector of the average firing rate at leaving
%           f_entering_all = vector of the average firing rate at entering
%           kinematics_leaving = vector of the average kinematics at leaving
%           kinematics_entering = vector of the average kinematics at entering
%           kinematics_leaving_std = vector of the standard dev in 
%               kinematics at leaving
%           kinematics_entering_std = vector of the standard dev in 
%               kinematics at entering
%   

df = f_orco.calcDeltaFR;
smoothWindow = 6;% smooth by 200 ms

df_leaving_aligned = cell(f_orco.nFly,1);
df_entering_aligned = cell(f_orco.nFly,1);
f_leaving_aligned = cell(f_orco.nFly,1);
f_entering_aligned = cell(f_orco.nFly,1);
kinematics_leaving_aligned = cell(f_orco.nFly,1);
kinematics_entering_aligned = cell(f_orco.nFly,1);
spd = f_orco.calcSpd;
for fly = 1:f_orco.nFly
    if strcmpi(type,'curv')
        kinematics = abs(f_orco.curv(fly,:)).*180./pi.*f_orco.fs;
    elseif strcmpi(type,'spd')
        kinematics = abs(spd(fly,:));
    end
    curr_f = f_orco.spk(fly,:);
    
    %leaving
    curr_df = -df(fly,:);
    [~,locs] = findpeaks(curr_df,'MinPeakHeight',thresh);
    locs(locs>(f_orco.nPt-delay)) = [];
    df_leaving_aligned{fly} = -curr_df(locs'+[-delay:delay]);
    f_leaving_aligned{fly} = curr_f(locs'+[-delay:delay]);
    kinematics_leaving_aligned{fly} = kinematics(locs'+[-delay:delay]);
    
    % entering
    curr_df = df(fly,:);
    [~,locs] = findpeaks(curr_df,'MinPeakHeight',thresh);
    locs(locs>(f_orco.nPt-delay)) = [];
    df_entering_aligned{fly} = curr_df(locs'+[-delay:delay]);
    f_entering_aligned{fly} = curr_f(locs'+[-delay:delay]);
    kinematics_entering_aligned{fly} = kinematics(locs'+[-delay:delay]);
    
end

kinematics_leaving_all = cell2mat(kinematics_leaving_aligned);
kinematics_entering_all = cell2mat(kinematics_entering_aligned);
df_leaving_all = nanmean(cell2mat(df_leaving_aligned));
df_entering_all = nanmean(cell2mat(df_entering_aligned));
f_leaving_all = nanmean(cell2mat(f_leaving_aligned));
f_entering_all = nanmean(cell2mat(f_entering_aligned));

% generate resampled data
nEmp = size(kinematics_leaving_all,1);
bootstrap_mean = zeros(100,delay*2+1);
for samp = 1:100
    X = datasample(1:nEmp,floor(nEmp./2),'Replace',false);
    bootstrap_sampled = kinematics_leaving_all(X,:);
    bootstrap_mean(samp,:) = smoothdata(nanmean(bootstrap_sampled),'movmean',smoothWindow);
end
kinematics_leaving_std = std(bootstrap_mean);

nEmp = size(kinematics_entering_all,1);
bootstrap_mean = zeros(100,delay*2+1);
for samp = 1:100
    X = datasample(1:nEmp,floor(nEmp./2),'Replace',false);
    bootstrap_sampled = kinematics_entering_all(X,:);
    bootstrap_mean(samp,:) = smoothdata(nanmean(bootstrap_sampled),'movmean',smoothWindow);
end
kinematics_entering_std = std(bootstrap_mean);

kinematics_leaving = smoothdata(nanmean(kinematics_leaving_all),'movmean',smoothWindow);
kinematics_entering = smoothdata(nanmean(kinematics_entering_all),'movmean',smoothWindow);

figure;set(gcf,'Position',[2 42 838 924]);
subplot(3,2,1);plot((-delay:delay)./30,kinematics_leaving);
xlabel('delay');ylabel(yl);title('Leaving')
subplot(3,2,2);plot((-delay:delay)./30,kinematics_entering);
xlabel('delay');ylabel(yl);title('Entering')
subplot(3,2,3);plot((-delay:delay)./30,df_leaving_all);
xlabel('delay');ylabel('spk/s^2');title('Leaving')
subplot(3,2,4);plot((-delay:delay)./30,df_entering_all);
xlabel('delay');ylabel('spk/s^2');title('Entering')
subplot(3,2,5);plot((-delay:delay)./30,f_leaving_all);
xlabel('delay');ylabel('spk/s');title('Leaving')
subplot(3,2,6);plot((-delay:delay)./30,f_entering_all);
xlabel('delay');ylabel('spk/s');title('Entering')

if singleTrial
    df_leaving_all = cellfun(@(x) mean(x),df_leaving_aligned,'UniformOutput',false);
    df_entering_all = cellfun(@(x) mean(x),df_entering_aligned,'UniformOutput',false);
    f_leaving_all = cellfun(@(x) mean(x),f_leaving_aligned,'UniformOutput',false);
    f_entering_all = cellfun(@(x) mean(x),f_entering_aligned,'UniformOutput',false);
    kinematics_leaving = cellfun(@(x) smoothdata(mean(x),'movmean',smoothWindow),kinematics_leaving_aligned,'UniformOutput',false);
    kinematics_entering = cellfun(@(x) smoothdata(mean(x),'movmean',smoothWindow),kinematics_entering_aligned,'UniformOutput',false);
    kinematics_leaving_std = cellfun(@(x) smoothdata(std(x)./sqrt(size(x,1)),'movmean',smoothWindow),kinematics_leaving_aligned,'UniformOutput',false);
    kinematics_entering_std = cellfun(@(x) smoothdata(std(x)./sqrt(size(x,1)),'movmean',smoothWindow),kinematics_entering_aligned,'UniformOutput',false);
% 
%     kinematics_leaving = smoothdata((kinematics_leaving_all),2,'movmean',smoothWindow);
%     kinematics_entering = smoothdata((kinematics_entering_all),2,'movmean',smoothWindow);
end

end