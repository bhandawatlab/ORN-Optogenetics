function [fNum] = plotKNNFormulation(f_orco,fNum)

xGrid = -150:15:150;
yGrid = 0:1:55;

% use curved walk as an example
param = 2;
m = f_orco.model.params{param};
allSpk = m.allSpk;
alldSpk = m.alldSpk;
allDat = m.allDat;

% get the KNN parameters
ratio = m.KNN.ratio(1:2);% not counting time
K = m.KNN.K;
Thresh = m.KNN.Thresh;

% KNN search 
X = [alldSpk, allSpk]./ratio;
[XX,YY] = ndgrid(xGrid,yGrid);
Y = [XX(:),YY(:)]./ratio;
[idx, D] = knnsearch(X,Y,'K',K);
mask = D<Thresh;

maskedDat = allDat(idx);
maskedDat(~mask) = nan;

%% Plotting functions
% plot the scatter plot of all tracks in 2D KNN space.
cNdx = {'g','r'};
figure(fNum);set(gcf,'Position',[2 42 838 924])
subplot(3,2,[1 2 3 4])
scatter(alldSpk,allSpk,25,'MarkerEdgeColor',[0, 0.4470, 0.7410],'linewidth',1);hold on
grid2cons = [0 20; 60 5];
tmpDat = cell(size(grid2cons,1),1);
for g = 1:size(grid2cons,1)
    YNdx = all(Y==grid2cons(g,:)./ratio,2);
    currPtNdx = idx(YNdx,:);
    currPtNdx = currPtNdx(D(YNdx,:)<Thresh);
    scatter(alldSpk(currPtNdx),allSpk(currPtNdx),25,cNdx{g},'filled')
    ellipse(Thresh.*ratio(1),Thresh.*ratio(2),0,grid2cons(g,1),grid2cons(g,2),'k')
    tmpDat{g} = maskedDat(YNdx,:);
end
plot([-150 150],[f_orco.spk(1) f_orco.spk(1)],'k')
plot([0 0],[0 45],'k')
legend({'all data','g1','g1 bound','g2', 'g2 bound'})
xlim([-150 150]);ylim([0 45]);
grid on
xlabel('df (spk/s^2)');ylabel('f (spk/s)')
title([f_orco.id ', ' m.state.state ', ' m.state.kin ', K=' num2str(K) ', Thresh=' num2str(Thresh)])
% plot sample distributions
for g = 1:2
    pHat = lognfit(tmpDat{g}(~isnan(tmpDat{g}))+eps);
    test_pdf = pdf('Lognormal',[0:0.5:30],pHat(1),pHat(2));
    
    subplot(3,2,4+g)
    histogram(tmpDat{g},[0:1:30],'Normalization','pdf','FaceColor','k');hold on;
    plot([0:0.5:30],test_pdf,'r','Linewidth',1);
    ylim([0 0.3]);
    xlabel('curved walk speed (mm/s)');
    ylabel('PDF')
    title(num2str(grid2cons(g,:)))
end
fNum = fNum+1;

% plot 3D representation
figure(fNum);set(gcf,'Position',[2 42 838 924]);hold on;view([45 45])
grid2cons = [0 20 30];
ratio3 = m.KNN.ratio;
ratio3(3) = ratio3(3)./f_orco.fs;% convert time ratio to seconds
for g = 1:size(grid2cons,1)
    [X,Y,Z] = ellipsoid(grid2cons(g,1),grid2cons(g,3),grid2cons(g,2),Thresh.*ratio3(1),Thresh.*ratio3(3),Thresh.*ratio3(2));
    h = surf(X,Y,Z);
    set(h,'FaceColor',cNdx{g}, ...
      'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none');hold on;
  l1 = [grid2cons(g,1) grid2cons(g,1)+Thresh.*ratio3(1)];
  l2 = [grid2cons(g,2) grid2cons(g,2)+Thresh.*ratio3(2)];
  l3 = [grid2cons(g,3) grid2cons(g,3)+Thresh.*ratio3(3)];
  plot3(l1, [grid2cons(g,3) grid2cons(g,3)], [grid2cons(g,2) grid2cons(g,2)],'k','Linewidth',2)
  plot3([grid2cons(g,1) grid2cons(g,1)], l3, [grid2cons(g,2) grid2cons(g,2)],'k','Linewidth',2)
  plot3([grid2cons(g,1) grid2cons(g,1)], [grid2cons(g,3) grid2cons(g,3)], l2,'k','Linewidth',2)
end
xlim([-150 150]);zlim([0 45]);ylim([0 90])
xlabel('df (spk/s^2)');zlabel('f (spk/s)');ylabel('time (s)')
view([45 45]);grid on
hold off;
fNum = fNum+1;

end


% figure;set(gcf,'Position',[2 42 838 924])
% subplot(3,2,[1 2 3 4])
% scatter(alldSpk,allSpk,25,'MarkerEdgeColor',[0, 0, 0],'linewidth',1);hold on
% grid2cons = [0, 20, 0; -90 10, 0; -30, 10, 0; 0, 0, 0; 60, 22, 0; 120, 15, 0];
% grid2cons = grid2cons(:,1:2);%no time;
% tmpDat = cell(size(grid2cons,1),1);
% %cNdx = {'g','r','b'};
% cNdx = varycolor(8);
% for g = 1:size(grid2cons,1)
%     YNdx = all(Y==grid2cons(g,:)./ratio,2);
%     currPtNdx = idx(YNdx,mask(YNdx,:));
%     scatter(alldSpk(currPtNdx),allSpk(currPtNdx),25,cNdx(g,:),'filled')
%     %scatter(dSpk_tmp(currPtNdx),spk_tmp(currPtNdx),25,cNdx(g,:),'filled');%cNdx{g}
%     %scatter(alldSpk(currPtNdx),allSpk(currPtNdx),cNdx{g})
%     %ellipse(Thresh./ratio(2),Thresh./ratio(1),0,grid2cons(g,1),grid2cons(g,2),'k')
%     ellipse(Thresh.*ratio(1),Thresh.*ratio(2),0,grid2cons(g,1),grid2cons(g,2),'k')
%     tmpDat{g} = maskedDat(YNdx,:);
%     tmpDat{g}(isnan(tmpDat{g})) = [];
% end
% plot([-150 150],[f_orco.spk(1) f_orco.spk(1)],'k')
% plot([0 0],[0 45],'k')
% legend({'all data','g1','g1 bound'})
% xlim([-150 150]);ylim([0 45]);
% grid on
% xlabel('df (spk/s^2)');ylabel('f (spk/s)')
% title([f_orco.id ', ' m.state.state ', ' m.state.kin ', K=' num2str(K) ', Thresh=' num2str(Thresh)])
