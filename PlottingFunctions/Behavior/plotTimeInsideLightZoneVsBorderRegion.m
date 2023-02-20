function [] = plotTimeInsideLightZoneVsBorderRegion(genAll,meta)
close all
%radialProb = cell(numel(genAll),1);
percentTimeSpentInA_before = cell(numel(genAll),1);
percentTimeSpentInA_during = cell(numel(genAll),1);

ttInA = @(x) sum(x<=1);
ttInB = @(x) sum(x>1 & x<=1.3);
for gen = 1:numel(genAll)
    
    load(strcat(string(meta.foldDataModel),'\',genAll{gen},'_',meta.d,meta.ext,'.mat'),'f_orco')
    stopSpd = 0;
    radialProb = f_orco.getRadialProb('H',meta.border,stopSpd,0);
    
    for fly = 1:f_orco.nFly
        beforeRPos = radialProb.before.yRaw{fly}.*f_orco.rBound;
        duringRPos = radialProb.during.yRaw{fly}.*f_orco.rBound;
        
        timeSpentInA_before = ttInA(beforeRPos);
        timeSpentInB_before = ttInB(beforeRPos);
        timeSpentInA_during = ttInA(duringRPos);
        timeSpentInB_during = ttInB(duringRPos);
        
        percentTimeSpentInA_before{gen}(fly,1) = timeSpentInA_before./(timeSpentInA_before+timeSpentInB_before);
        percentTimeSpentInA_during{gen}(fly,1) = timeSpentInA_during./(timeSpentInA_during+timeSpentInB_during);
    end
    nFly(gen) = f_orco.nFly;
end

% define basic paramters
lims = [0 1];                      % y limits to set (example: [0 1])
isPaired = 'N';                 % equivalent to paired t-test
circleSize = 60;                % size of scatter plot points
barstate = 'off';                % 'On' = bar graph, 'Off' = scatter plot
subplots = [3,2,1,1]; % ysbplt, xsbplt, sbpltN, fig
xSubplots = repmat(100./subplots(1),1,subplots(1));
ySubplots = repmat(100./subplots(2),1,subplots(2));

figure(subplots(4));set(gcf,'Position',[842 42 838 924])
p = panel();
p.pack(xSubplots, ySubplots);


for gen = 1:numel(genAll)
    for gen2 = (gen+1):numel(genAll)
        if subplots(3)>6
            subplots(4) = subplots(4)+1;
            figure(subplots(4));set(gcf,'Position',[842 42 838 924])
            subplots(3) = 1;
        end
        
        data = 1-[percentTimeSpentInA_during{gen};percentTimeSpentInA_during{gen2}];% 1-% time spent in A
        g = [repmat(genAll(gen),nFly(gen),1);repmat(genAll(gen2),nFly(gen2),1)];
        [ss,p,~] = dabest3(data+eps*rand(size(data)),g,p,[],lims,isPaired,circleSize,barstate,subplots);
        labelAxis(p,subplots,'Percent time','1<r<1.3cm/(r<1.3cm)');
        subplots(3) = subplots(3)+1;
    end
    sgtitle('Percent time spent in border (1-1.3cm radial dist/(<1.3 cm radial dist))')
end
printFigures(1:subplots(4),[char(meta.plotFold) '\GeneralBehavior\'],'TimeInLightBorderVsBorderAndInnerLightZone');

% figure;set(gcf,'Position',[2 42 838 924]);
% k = 1;
% opts.plotData = false;
% opts.yl = [0 1];
% for gen = 1:numel(genAll)
%     for gen2 = (gen+1):numel(genAll)
%         data = nan(max([nFly(gen) nFly(gen2)]),2);
%         data(1:nFly(gen),1) = percentTimeSpentInA_during{gen};
%         data(1:nFly(gen2),2) = percentTimeSpentInA_during{gen2};
%         opts.xLabels = genAll([gen gen2]);
%         subplot(3,2,k)
%         violinPlotsStats(data,opts)
%         ylabel('Percent time');
%         xtickangle(10)
%         k = k+1;
%     end
%     sgtitle('Percent time spent in region A (<1cm radial dist/(1-1.3 cm radial dist))')
% end

end

function [] = labelAxis(p,subplots,label,titl)
[J,I] = ind2sub(subplots(2:-1:1),subplots(3));
try
    p(I,J).select();
    ylabel(label);
    title(titl,'Interpreter','none')
    xtickangle(10)
catch
    p(I,J,1,1).select();
    ylabel(label);
    title(titl,'Interpreter','none')
    p(I,J,2,1).select();
    ylabel(['delta ' label]);
    xtickangle(10)
end

end

function [] = printFigures(fNum,folder,figTitle)
if ~exist(folder, 'dir')
    mkdir(folder)
end

if numel(fNum)==1
    fRange = 1:max(fNum-1,1);
else
    fRange = fNum;
end
fRange(isnan(fRange)) = [];
if exist([folder figTitle '.ps'], 'file')==2
  delete([folder figTitle '.ps']);
end

for f = fRange
    figure(f);
    print('-bestfit','-painters','-dpsc2',[folder figTitle '.ps'],'-loose','-append');
end
ps2pdf('psfile', [folder figTitle '.ps'], 'pdffile', [folder figTitle '.pdf'], 'gspapersize', 'letter',...
    'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
    'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
    'gslibpath','C:\Program Files\gs\gs9.50\lib');
end