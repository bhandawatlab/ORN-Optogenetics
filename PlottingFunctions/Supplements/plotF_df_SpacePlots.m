function [fNum] = plotF_df_SpacePlots(f_orcoAll,fNum,border,Intensity2VoltageFile)
alldSpk = [];allSpk = [];allRPos = [];
alldSpkSmooth = [];allSpkSmooth = [];
dF_thresh = 20;F_thresh = 15;

for g = 1:numel(f_orcoAll)
    f_orco = f_orcoAll{g};
    gen = f_orco.id;
    hist = ceil(0.2.*f_orco.fs)-1;%200 ms history
    
    % calculate the different states
    tic;[cross] = f_orco.getCrossingBaseline('H',border,4.8,0);toc;
    dSpk = f_orco.calcDeltaFR;
    baselineSpk = f_orco.spk(1);
    
    dFSmooth = dSpk(:,1:end-hist);
    FSmooth = f_orco.spk(:,1:end-hist);
    for i = 2:hist+1
        dFSmooth = dSpk(:,i:end-(hist+1)+i)+dFSmooth;
        FSmooth = f_orco.spk(:,i:end-(hist+1)+i)+FSmooth;
    end
    dFSmooth = [zeros(size(dFSmooth,1),hist),dFSmooth./(hist+1)];% average
    FSmooth = [f_orco.spk(1).*ones(size(FSmooth,1),hist),FSmooth./(hist+1)];% average
    
    
    if g == 1
        flySamp = [7 2];iSamp = [10 3];

        figure(fNum);set(gcf,'Position',[2 42 838 924])
        for idx = 1:numel(flySamp)
            fly = flySamp(idx);

            ndx = cross.track{fly,iSamp(idx)}(4,:);
            n = ndx(end);
            while abs(dSpk(fly,n))>0.001 || f_orco.spk(fly,n)<4.7
                n = n+1;
            end
            ndx = ndx(1)-4:n;

            tt = (1:1:numel(ndx))./f_orco.fs;
            nearestInt = round2NearestInterval(tt(end),1);

            %plot(dSpk(fly,ndx), f_orco.spk(fly,ndx));
            %plot3(dSpk(fly,ndx), f_orco.spk(fly,ndx),(1:1:numel(ndx))./f_orco.fs);
            subplot(3,2,idx);
            plot(dSpk(fly,ndx), f_orco.spk(fly,ndx),'k','LineWidth',2);hold on;
            plot(dFSmooth(fly,ndx), FSmooth(fly,ndx),'Color',[0.5 0.5 0.5],'LineWidth',2);
            plot([-150 150],[baselineSpk baselineSpk],'k','LineWidth',1)
            plot([0 0],[0 45],'k','LineWidth',1);
            plot([dF_thresh dF_thresh],[0 45],'r');
            plot(-[dF_thresh dF_thresh],[0 45],'r')
            xlim([-150 150]);ylim([0 45])
            xlabel('delta spike rate (spk/s^2)');ylabel('spike rate (spk/s)')
            subplot(3,2,idx+2);
            plot3(tt, dSpk(fly,ndx), f_orco.spk(fly,ndx),'k','LineWidth',2);hold on;
            plot3(tt, dFSmooth(fly,ndx), FSmooth(fly,ndx),'Color',[0.5 0.5 0.5],'LineWidth',2);
            plot3([0 nearestInt],[0 0],[baselineSpk baselineSpk],'k','LineWidth',1)
            xlim([0 nearestInt]);ylim([-150 150]);zlim([0 45])
            xlabel('time (s)');ylabel('delta spike rate (spk/s^2)');zlabel('spike rate (spk/s)')
            subplot(3,2,idx+4);yyaxis left;plot(tt,f_orco.spk(fly,ndx),'LineWidth',2);hold on;
            plot(tt,FSmooth(fly,ndx),'Color',[0.5 0.5 0.5],'LineWidth',2)
            xlim([0 nearestInt]);ylim([0 45]);ylabel('spike rate')
            yyaxis right;plot(tt,f_orco.rH(fly,ndx),'LineWidth',2)
            xlim([0 nearestInt]);ylim([0 4]);ylabel('radial position')
            xlabel('time');
        end
        suptitle(['Sample Trajectory in F/dF space: ' gen])
        fNum = fNum+1;
    end
    try
        alldSpkSmooth = [alldSpkSmooth;dFSmooth];%[alldSpk;dSpk];
        allSpkSmooth = [allSpkSmooth;FSmooth];%[allSpk;f_orco.spk];
        alldSpk = [alldSpk;dSpk];
        allSpk = [allSpk;f_orco.spk];
        allRPos = [allRPos;f_orco.rH];
    catch
        alldSpkSmooth = [alldSpkSmooth;[zeros(size(dFSmooth,1),1),dFSmooth]];%[alldSpk;dSpk];
        allSpkSmooth = [allSpkSmooth;[zeros(size(dFSmooth,1),1),FSmooth]];%[allSpk;f_orco.spk];
        alldSpk = [alldSpk;[zeros(size(dFSmooth,1),1),dSpk]];
        allSpk = [allSpk;[zeros(size(dFSmooth,1),1),f_orco.spk]];
        allRPos = [allRPos;[zeros(size(dFSmooth,1),1),f_orco.rH]];
    end
end

badNdx = alldSpk==0 & allSpk==baselineSpk;
alldSpkSmooth(badNdx) = [];%[alldSpk;dSpk];
allSpkSmooth(badNdx) = [];%[allSpk;f_orco.spk];
alldSpk(badNdx) = [];
allSpk(badNdx) = [];

ndx = datasample(1:numel(alldSpk),10000);

figure(fNum);set(gcf,'Position',[2 42 838 924])
subplot(2,1,1);
scatter(alldSpk(ndx),abs(allSpk(ndx)),[],[0.5 0.5 0.5]);hold on;
plot([-20 20],[F_thresh F_thresh],'k','LineWidth',1)
plot([-20 20],[baselineSpk baselineSpk],'k','LineWidth',1)
plot([0 0],[0 50],'k','LineWidth',1);
plot([dF_thresh dF_thresh],[0 45],'r');
plot(-[dF_thresh dF_thresh],[0 45],'r')
xlim([-400 400])
subplot(2,1,2);
scatter(alldSpkSmooth(ndx),abs(allSpkSmooth(ndx)),[],[0.5 0.5 0.5]);hold on;
plot([-20 20],[F_thresh F_thresh],'k','LineWidth',1)
plot([-20 20],[baselineSpk baselineSpk],'k','LineWidth',1)
plot([0 0],[0 50],'k','LineWidth',1);
plot([dF_thresh dF_thresh],[0 45],'r');
plot(-[dF_thresh dF_thresh],[0 45],'r');
xlim([-250 250])
fNum = fNum+1;

% h = histcounts(abs(alldSpk(abs(alldSpk)>0)),[0:1:200]);
% figure(fNum);set(gcf,'Position',[2 42 838 924])
% subplot(3,1,1);
% plot([0:1:200],[0 cumsum(h)./sum(h)]);hold on;plot([dF_thresh,dF_thresh],[0 max(1)]);
% xlabel('|delta firing rate (spk/s^2)|');ylabel('CDF');
% N = histcounts2(alldSpk(abs(alldSpk)>dF_thresh),allSpk(abs(alldSpk)>dF_thresh),[-200:10:200],[0:2:50]);
% subplot(3,1,[2 3]);imagesc([-200:10:200],[0:2:50],N');colormap(hot);set(gca,'YDir','normal');colorbar
% xlabel('delta firing rate (spk/s^2)');ylabel('firing rate (spk/s)');
% fNum = fNum+1;
% 
% 
% 
% if ~isempty(Intensity2VoltageFile)
%     I = f_orco.model.IntensitySpace.I;
%     x = f_orco.model.IntensitySpace.x;
%     %convIV = f_orco.model.IntensitySpace.convIV;
%     load(Intensity2VoltageFile,'p','convIV')
%     rPos = allRPos;
%     for i = 1:size(allRPos,2)
%         [V_out(:,i),I_out(:,i)] = location2Intensity(rPos(:,i),x,I,convIV);
%     end
%     dI_out = gradient(I_out).*f_orco.fs;
% 
%     h = histcounts(abs(dI_out(abs(dI_out)>0)),[0:0.2:50]);
%     figure(fNum);set(gcf,'Position',[2 42 838 924]);subplot(3,1,1);
%     plot([0:0.2:50],[0 cumsum(h)./sum(h)]);hold on;plot([3,3],[0 max(1)]);
%     xlabel('|delta Intensity (1/s)|');ylabel('CDF');
%     N = histcounts2(dI_out(abs(dI_out)>3),I_out(abs(dI_out)>3),[-50:1:50],[0:0.15:3.5]);
%     subplot(3,1,[2 3]);imagesc([-50:0.2:50],[0:0.1:3.5],N');colormap(hot);set(gca,'YDir','normal');colorbar
%     xlabel('delta Intensity (1/s)');ylabel('Intensity (mW/cm^2)');
%     fNum = fNum+1;
% end

end
