function [df_leaving_all,df_entering_all,f_leaving_all,f_entering_all,...
    turnProb_leaving,turnProb_entering,turnProb_leaving_std,turnProb_entering_std] = ...
    getDf_aligned_turnRates(f_orco,thresh,delay,turnNdxAbs,transitionsOnly,yl,fs)
%getDf_aligned_turnRates  calculates the firing rate, change in firing
%   rate, and the turn probability when leaving or entering (aligned by 
%   peak in change in firing rate)
%
%   Inputs: f_orco = fly object
%           thresh = threshold for peak (or trough) in delta firing rate
%           delay = time window to consider before and after the peak
%           turnNdxAbs = fly x 1 cell array where each cell indicates the
%               times at which a sharp turn transition occurs for the fly
%           transitionsOnly = true/false (true = transition prob, false =
%               prob of being in sharp turn)
%           yl = ylabel (for plotting purposes)
%           fs = sampling rate
%
%   Output: df_leaving_all = vector of the average change in firing rate
%               at leaving (-delay to +delay)
%           df_entering_all = vector of the average change in firing rate
%               at entering
%           f_leaving_all = vector of the average firing rate at leaving
%           f_entering_all = vector of the average firing rate at entering
%           turnProb_leaving = vector of the turn prob at leaving
%           turnProb_entering = vector of the turn prob at entering
%           turnProb_leaving_std = vector of the standard dev in turn prob
%               at leaving
%           turnProb_entering_std = vector of the standard dev in turn prob
%               at entering
%   

df = f_orco.calcDeltaFR;
stNdx = find(strcmpi(f_orco.states.key  ,'sharp turns'));

df_leaving_aligned = cell(f_orco.nFly,1);
df_entering_aligned = cell(f_orco.nFly,1);
f_leaving_aligned = cell(f_orco.nFly,1);
f_entering_aligned = cell(f_orco.nFly,1);
turn_leaving_aligned = cell(f_orco.nFly,1);
turn_entering_aligned = cell(f_orco.nFly,1);
for fly = 1:f_orco.nFly
    if transitionsOnly
        % ---------- transition into ----------
        turn_binary = zeros(1,f_orco.nPt);
        turn_binary(turnNdxAbs{fly}) = 1;
    else
        % ---------- in vs out of a sharp turn ----------
        turn_binary = f_orco.states.ndx(fly,:)==stNdx;
    end
    curr_f = f_orco.spk(fly,:);
    
    %leaving
    curr_df = -df(fly,:);
    [~,locs] = findpeaks(curr_df,'MinPeakHeight',thresh);
    locs(locs>(f_orco.nPt-delay)) = [];
    df_leaving_aligned{fly} = -curr_df(locs'+[-delay:delay]);
    f_leaving_aligned{fly} = curr_f(locs'+[-delay:delay]);
    turn_leaving_aligned{fly} = turn_binary(locs'+[-delay:delay]);
    
    % entering
    curr_df = df(fly,:);
    [~,locs] = findpeaks(curr_df,'MinPeakHeight',thresh);
    locs(locs>(f_orco.nPt-delay)) = [];
    df_entering_aligned{fly} = curr_df(locs'+[-delay:delay]);
    f_entering_aligned{fly} = curr_f(locs'+[-delay:delay]);
    turn_entering_aligned{fly} = turn_binary(locs'+[-delay:delay]);
    
end

tmp_df_leaving = cell2mat(df_leaving_aligned);
tmp_df_entering = cell2mat(df_entering_aligned);
tmp_f_leaving = cell2mat(f_leaving_aligned);
tmp_f_entering = cell2mat(f_entering_aligned);


turn_leaving_all = cell2mat(turn_leaving_aligned);
turn_entering_all = cell2mat(turn_entering_aligned);
df_leaving_all = nanmean(tmp_df_leaving);
df_entering_all = nanmean(tmp_df_entering);
f_leaving_all = nanmean(tmp_f_leaving);
f_entering_all = nanmean(tmp_f_entering);


% df_leaving_sem = nanstd(tmp_df_leaving)./size(tmp_df_leaving,1);
% df_entering_sem = nanstd(tmp_df_entering)./size(tmp_df_entering,1);
% f_leaving_sem = nanstd(tmp_f_leaving)./size(tmp_f_leaving,1);
% f_entering_sem = nanstd(tmp_f_entering)./size(tmp_f_entering,1);

df_leaving_std = nanstd(tmp_df_leaving);
df_entering_std = nanstd(tmp_df_entering);
f_leaving_std = nanstd(tmp_f_leaving);
f_entering_std = nanstd(tmp_f_entering);

% generate resampled data
nEmp = size(turn_leaving_all,1);
bootstrap_mean = zeros(100,delay*2+1);
for samp = 1:100
    X = datasample(1:nEmp,floor(nEmp./2),'Replace',false);
    bootstrap_sampled = turn_leaving_all(X,:);
    bootstrap_mean(samp,:) = smoothdata(nanmean(bootstrap_sampled),'movmean',30);
end
turnProb_leaving_std = std(bootstrap_mean);

nEmp = size(turn_entering_all,1);
bootstrap_mean = zeros(100,delay*2+1);
for samp = 1:100
    X = datasample(1:nEmp,floor(nEmp./2),'Replace',false);
    bootstrap_sampled = turn_entering_all(X,:);
    bootstrap_mean(samp,:) = smoothdata(nanmean(bootstrap_sampled),'movmean',30);
end
turnProb_entering_std = std(bootstrap_mean);

turnProb_leaving = smoothdata(nanmean(turn_leaving_all),'movmean',30);
turnProb_entering = smoothdata(nanmean(turn_entering_all),'movmean',30);

figure;set(gcf,'Position',[2 42 838 924]);
subplot(3,2,1);shadedErrorBar((-delay:delay)./fs,turnProb_leaving,turnProb_leaving_std);
xlabel('delay');ylabel(yl);title('Leaving')
subplot(3,2,2);shadedErrorBar((-delay:delay)./fs,turnProb_entering,turnProb_entering_std);
xlabel('delay');ylabel(yl);title('Entering')
subplot(3,2,3);shadedErrorBar((-delay:delay)./fs,df_leaving_all,df_leaving_std);
xlabel('delay');ylabel('spk/s^2');title('Leaving')
subplot(3,2,4);shadedErrorBar((-delay:delay)./fs,df_entering_all,df_entering_std);
xlabel('delay');ylabel('spk/s^2');title('Entering')
subplot(3,2,5);shadedErrorBar((-delay:delay)./fs,f_leaving_all,f_leaving_std);
xlabel('delay');ylabel('spk/s');title('Leaving')
subplot(3,2,6);shadedErrorBar((-delay:delay)./fs,f_entering_all,f_entering_std);
xlabel('delay');ylabel('spk/s');title('Entering')

% figure;set(gcf,'Position',[2 42 838 924]);
% subplot(3,2,1);plot((-delay:delay)./fs,turnProb_leaving);
% xlabel('delay');ylabel(yl);title('Leaving')
% subplot(3,2,2);plot((-delay:delay)./fs,turnProb_entering);
% xlabel('delay');ylabel(yl);title('Entering')
% subplot(3,2,3);plot((-delay:delay)./fs,df_leaving_all);
% xlabel('delay');ylabel('spk/s^2');title('Leaving')
% subplot(3,2,4);plot((-delay:delay)./fs,df_entering_all);
% xlabel('delay');ylabel('spk/s^2');title('Entering')
% subplot(3,2,5);plot((-delay:delay)./fs,f_leaving_all);
% xlabel('delay');ylabel('spk/s');title('Leaving')
% subplot(3,2,6);plot((-delay:delay)./fs,f_entering_all);
% xlabel('delay');ylabel('spk/s');title('Entering')

end