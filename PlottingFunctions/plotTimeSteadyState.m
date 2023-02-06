function [fNum] = plotTimeSteadyState(f_orcoAll,border,fNum)
% time spent at each steady state and for entry/exit
dur_baseline = cell(4,numel(f_orcoAll));
for i = 1:numel(f_orcoAll)
    fs = f_orcoAll{i}.fs;
    spk = f_orcoAll{i}.spk;
    dSpk = f_orcoAll{i}.calcDeltaFR;
    baseline = mean(spk(:,30*fs:60*fs),'all');
    fe = f_orcoAll{i}.getFirstEntry('H',border);
    
    % separate out into each time point into before first entry (FE), below
    % baseline firing rate, baseline firing rate after FE, and above
    % baseline firing rate
    % enter/leaving
    condNdx = 4.*ones(size(spk));
    % baseline SS
    condNdx(abs((spk-baseline))<1 & abs(dSpk)<3) = 1;
    %condNdx(abs((spk-baseline))<0.001 & dSpk==0) = 1;
    % inhibition SS
    condNdx(spk<(baseline-1) & abs(dSpk)<3) = 2;
    % inside ss
    condNdx(spk>(baseline+1) & abs(dSpk)<3) = 3;
%     % enter/leaving
%     condNdx(spk>(baseline+0.001) & abs(dSpk)<0.1) = 4;
    dur_baseline{4,i} = [];
    
    for j = 1:f_orcoAll{i}.nFly
        for cond = 1:4
            condNdx(j,1:fe(j)-1) = 0;
            [startNdx,endNdx,type] = startEndSeq(condNdx(j,:)==cond);
            startNdx = startNdx(type);
            endNdx = endNdx(type);
            
            dur_baseline{cond,i} = [dur_baseline{cond,i},endNdx-startNdx+1];
        end
    end
    fNames{i} = erase(f_orcoAll{i}.id,' Retinal');
end

noiseThresh = 30;
for cond = 1:4
    for i = 1:numel(f_orcoAll)
        tmp = dur_baseline{cond,i};
        tmp_noNoise = tmp(tmp>noiseThresh);
        totTime_cond(cond,i) = sum(tmp);
        totTime_cond_nonoise(cond,i) = sum(tmp_noNoise);
    end
end
X = categorical(fNames);

% percent "noise"
totTime_cond_noise = (totTime_cond-totTime_cond_nonoise);
totTime_noise_baseline_enterLeave = [sum(totTime_cond_noise);...
    sum(totTime_cond_nonoise(1:3,:));totTime_cond_nonoise(4,:)];

figure(fNum);
set(gcf,'Position',[438   119   559   783]);
subplot(7,1,[1:2]);
barh(X,(totTime_cond./sum(totTime_cond))','stacked');
xlim([0 1])
legend({'baseline SS','inhibition SS','inside ss','enter/leave'})
title(['All tracks'])
subplot(7,1,[3:4]);
barh(X,(totTime_cond_nonoise./sum(totTime_cond_nonoise))','stacked')
xlim([0 1])
legend({'baseline SS','inhibition SS','inside ss','enter/leave'})
title(['Only tracks> ' num2str(noiseThresh./30) ' s'])
subplot(7,1,[5:6]);
barh(X,(totTime_noise_baseline_enterLeave./sum(totTime_noise_baseline_enterLeave))','stacked')
xlim([0 1])
legend({['<' num2str(num2str(noiseThresh./30)) 's tracks'],...
    'Steady State','enter/leave'})
title(['Percent of tracks enter/leave vs SS vs noise'])

condition = {'baseline SS';'inhibition SS';'inside SS';'enter/leave'};
%f = [1;-1;1;nan];
%df = [3;3;3;nan];

f = {'abs((f-baseline))<1';'f<(baseline-1)';'f>(baseline+1)';'otherwise'};
df = {'abs(df)<3';'abs(df)<3';'abs(df)<3';'otherwise'};
T = table(condition,f,df,'VariableNames',{'condition','f thresh','df thresh'});

% Convert Table to cell to char array
tableCell = [T.Properties.VariableNames; table2cell(T)]; 
tableCell(cellfun(@isnumeric,tableCell)) = cellfun(@num2str, tableCell(cellfun(@isnumeric,tableCell)),'UniformOutput',false); 
tableChar = splitapply(@strjoin,pad(tableCell),[1;2;3;4;5]);
% Add axes (not visible) & text (use a fixed width font)
axes('position',[.1,.1,.8,.2], 'Visible','off')
text(.1,.4,tableChar,'VerticalAlignment','Top','HorizontalAlignment','Left','FontName','Consolas');

end