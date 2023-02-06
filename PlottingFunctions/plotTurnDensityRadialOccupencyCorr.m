function [fNum] = plotTurnDensityRadialOccupencyCorr(radProbAllGen,turnDensAllGen,genAll,fNum)

posDiff = [];turnDiff = [];
for g = 1:numel(genAll)
    posDiff(g,:) = radProbAllGen{1,g}.during.y-radProbAllGen{2,g}.during.y;
    turnDiff(g,:) = turnDensAllGen{1,g}.Dur-turnDensAllGen{2,g}.Dur;
end

[rhoCorr,~] = corr(posDiff',turnDiff');
D = diag(rhoCorr);

posDiffPosOnly = posDiff;
posDiffPosOnly(posDiffPosOnly<0) = 0;
turnDiffPosOnly = turnDiff;
turnDiffPosOnly(turnDiffPosOnly<0) = 0;
turnDiffNegOnly = turnDiff;
turnDiffNegOnly(turnDiffNegOnly>0) = 0;


figure(fNum);set(gcf,'Position',[2 42 838 924])
xx = (turnDensAllGen{1,1}.x(1:end-1)+turnDensAllGen{1,1}.x(2:end))./2.*4;
for g = 1:numel(genAll)
    subplot(5,3,g);
    yyaxis left;plot(xx,posDiff(g,:));ylim([-0.2 0.2]);ylabel('radial loc diff');
    xlabel('radial position');xticks([0:0.5:4])
    yyaxis right;plot(xx,turnDiff(g,:));ylim([-0.1 0.1]);ylabel('Turn density diff')
    title(genAll{g})
end
fNum = fNum+1;

figure(fNum);set(gcf,'Position',[2 42 838 924])
xx = (turnDensAllGen{1,1}.x(1:end-1)+turnDensAllGen{1,1}.x(2:end))./2.*4;
for g = 1:numel(genAll)
    subplot(5,3,g);
    scatter(turnDiff(g,:),posDiff(g,:));
    axis square;
    xlim([-0.2 0.2]);ylim([-0.2 0.2])
    h = lsline;xlim([-0.2 0.2]);ylim([-0.2 0.2])
    h.Color = 'r';
    h.LineWidth = 1;
    ylabel('radial loc diff');
    xlabel('Turn density diff');%xticks([0:0.5:4])
    title([genAll{g} ' Rsq=' num2str(D(g).^2)])
end
fNum = fNum+1;

figure(fNum);set(gcf,'Position',[2 42 838 924])
hold on;
x2 = [xx, fliplr(xx)];
inBetween = [posDiff(1,:), fliplr(fliplr(zeros(size(xx))))];
fill(x2, inBetween, 'g');
plot(xx,posDiff(1,:),'k');ylim([-0.2 0.2]);
plot(xx,zeros(size(xx)),'k')
xticks([0:0.5:4])
xlabel('radial position');ylabel('Radial Occupancy diff');
title(genAll{1})
fNum = fNum+1;

figure(fNum);set(gcf,'Position',[2 42 838 924])
hold on;
x2 = [xx, fliplr(xx)];
inBetween = [turnDiff(1,:), fliplr(fliplr(zeros(size(xx))))];
fill(x2, inBetween, 'g');
plot(xx,turnDiff(1,:),'k');ylim([-0.1 0.1]);
plot(xx,zeros(size(xx)),'k')
xticks([0:0.5:4])
xlabel('radial position');ylabel('Turn density diff');
title(genAll{1})
fNum = fNum+1;

figure(fNum);set(gcf,'Position',[2 42 838 924])
md1 = fitlm(sum(posDiffPosOnly,2),D);
subplot(3,2,1);h = plot(md1);h(1).Marker = 'o';
xlim([0.1 0.6]);ylim([0 1]);title(['r-sq=' num2str(md1.Rsquared.Ordinary)])
xlabel('total positive pos diff');ylabel('Correlation')

md1 = fitlm(sum(turnDiffPosOnly,2),D);
subplot(3,2,2);h = plot(md1);h(1).Marker = 'o';
xlim([0.1 0.35]);ylim([0 1]);title(['r-sq=' num2str(md1.Rsquared.Ordinary)])
xlabel('total positive turn diff');ylabel('Correlation')

md1 = fitlm(max(posDiffPosOnly,[],2),D);
subplot(3,2,3);h = plot(md1);h(1).Marker = 'o';
xlim([0 0.2]);ylim([0 1]);title(['r-sq=' num2str(md1.Rsquared.Ordinary)])
xlabel('max positive pos diff');ylabel('Correlation')

md1 = fitlm(max(turnDiffPosOnly,[],2),D);
subplot(3,2,4);h = plot(md1);h(1).Marker = 'o';
xlim([0 0.15]);ylim([0 1]);title(['r-sq=' num2str(md1.Rsquared.Ordinary)])
xlabel('max positive turn diff');ylabel('Correlation')

pd1 = fitdist(D,'kernel');
y = 0:0.1:1;
xLeft = -pdf(pd1,y);
xRight = pdf(pd1,y);
x = pdf(pd1,D).*(rand(size(D))-0.5);

subplot(3,2,5);scatter(x,D);hold on;
plot(xLeft,y,'Color','r','LineStyle','-')
plot(xRight,y,'Color','r','LineStyle','-')
ylabel('correlation')
subplot(3,2,6);bar(0,mean(D),'w');hold on
er = errorbar(0,mean(D),std(D),std(D));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
scatter(x./2,D,'');
ylabel('correlation')
hold off
fNum = fNum+1;


end


