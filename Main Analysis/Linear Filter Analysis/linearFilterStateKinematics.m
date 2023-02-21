function [] = linearFilterStateKinematics(gen,meta)

close all
if ~exist(strcat(string(meta.plotFold),'/LinearFilterAnalysis/'), 'dir')
    mkdir(strcat(string(meta.plotFold),'/LinearFilterAnalysis/'));
end
psFileName = strcat(string(meta.plotFold),'/LinearFilterAnalysis/LinearFilterStateKinematics.ps');
pdfFileName = strcat(string(meta.plotFold),'/LinearFilterAnalysis/LinearFilterStateKinematics.pdf');
if exist(psFileName, 'file')==2
    delete(psFileName);
end

load(strcat(string(meta.foldDataModel),'\',gen,'_',meta.d,meta.ext,'.mat'),'f_orco')

ntfilt_all = round(5.*f_orco.fs);%seconds;

stateHist = cell(4,1);
b_all = cell(4,1);
rho_all = cell(4,1);
eta_all = cell(4,1);
for kin = 1:4
    
    stateFly = f_orco.model.params{1,kin}.allFly;
    stateTime = f_orco.model.params{1,kin}.allTime;
    stateKin = f_orco.model.params{1,kin}.allDat;
    
    stateHist{kin} = cell(f_orco.nFly,1);
    for fly = 1:f_orco.nFly
        currSpk = f_orco.spk(fly,:);
        ndx = stateTime(stateFly==fly)+(-ntfilt_all:1:-1)+1;
        stateHist{kin}{fly,1} = [ones(sum(stateFly==fly),1) currSpk(ndx)];
    end
    spkHist_all = cell2mat(stateHist{kin});
    
    tic;[U,s,V] = csvd(spkHist_all,[]);toc
    %tic;[U,s,V] = svd(XStimNew);%toc
    % solve for the ridge regression
    nLamb = 30;
    lamb = 2.^((1:nLamb)-15);
    b = zeros(ntfilt_all+1,nLamb);rho = zeros(1,nLamb);eta = zeros(1,nLamb);
    for i = 1:nLamb
        [b(:,i),rho(i),eta(i)] = tikhonov(U,s,V,stateKin,lamb(i),zeros(ntfilt_all+1,1));
    end
    b_all{kin} = b;
    rho_all{kin} = b;
    eta_all{kin} = b;
end

%%
%load('tmp_LinearFilterStateKinematics.mat')
close all
%

figure;set(gcf,'Position',[2 42 838 924]);
for kin = 1:4
    spkHist_all = cell2mat(stateHist{kin});
    closestFit = max(spkHist_all*b_all{kin}(:,1),0);
    subplot(2,2,kin);
    scatter((f_orco.model.params{1,kin}.allDat),closestFit,'k');hold on;
    roundMax = round2NearestInterval(max(f_orco.model.params{1,kin}.allDat),5);
    plot(([0 roundMax]),([0 roundMax]),'r');
    xlim([0 roundMax]);ylim([0 roundMax])
    ax = gca;
    set(ax,'YTick',get(ax,'XTick'));
    xlabel('empirical');ylabel('predicted (h*x)');axis tight
    title([f_orco.model.params{1,kin}.state.state f_orco.model.params{1,kin}.state.kin])
end

k = 64;
figure;set(gcf,'Position',[2 42 838 924]);
for kin = 1:4
    stateKin = log(f_orco.model.params{1,kin}.allDat);
    stateSpk = f_orco.model.params{1,kin}.allSpk;
    statedSpk = f_orco.model.params{1,kin}.alldSpk;
    
    [stateSpk_sort,ndx] = sort(stateSpk);
    stateKin_sort = stateKin(ndx);
    
    moveMean = smooth(stateSpk_sort, stateKin_sort,k);
    %moveMean2 = smooth(stateKin_sort,64);
    moveStd = movstd(stateKin_sort,k);
    
    subplot(2,2,kin);
    scatter(stateSpk, stateKin,'k');hold on;
    plot(stateSpk_sort,moveMean,'r-','LineWidth',2)
    plot(stateSpk_sort,moveMean+moveStd,'r--','LineWidth',2);
    plot(stateSpk_sort,moveMean-moveStd,'r--','LineWidth',2);
    xlabel('200 ms mean spk (h_p*x)');ylabel('true');
    title([f_orco.model.params{1,kin}.state.state f_orco.model.params{1,kin}.state.kin])
end

k = 64;
figure;set(gcf,'Position',[2 42 838 924]);
for kin = 1:4
    stateKin = log(f_orco.model.params{1,kin}.allDat);
    stateSpk = f_orco.model.params{1,kin}.allSpk;
    statedSpk = f_orco.model.params{1,kin}.alldSpk;
    
    [stateSpk_sort,ndx] = sort(statedSpk);
    stateKin_sort = stateKin(ndx);
    
    moveMean = smooth(stateSpk_sort, stateKin_sort,k);
    %moveMean2 = smooth(stateKin_sort,64);
    moveStd = movstd(stateKin_sort,k);
    
    subplot(2,2,kin);
    scatter(statedSpk, stateKin,'k');hold on;
    plot(stateSpk_sort,moveMean,'r-','LineWidth',2)
    plot(stateSpk_sort,moveMean+moveStd,'r--','LineWidth',2);
    plot(stateSpk_sort,moveMean-moveStd,'r--','LineWidth',2);
    xlabel('delta spk (g_p*x)');ylabel('true');
    title([f_orco.model.params{1,kin}.state.state f_orco.model.params{1,kin}.state.kin])
end

for f = 1:get(gcf,'Number')
    figure(f);
    print('-painters','-dpsc2',psFileName,'-loose','-append');
end
try
    ps2pdf('psfile',psFileName, 'pdffile', ...
        pdfFileName, 'gspapersize', 'letter',...
        'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
        'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
        'gslibpath','C:\Program Files\gs\gs9.50\lib');
catch
    disp('No ghostscript available. Please install ghostscript or ')
    disp('change path to ghostscript in line 120 of linearFilterStateKinematics.m')
end
end