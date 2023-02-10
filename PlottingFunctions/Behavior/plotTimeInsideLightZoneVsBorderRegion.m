function [] = plotTimeInsideLightZoneVsBorderRegion(genAll,meta)

%radialProb = cell(numel(genAll),1);
percentTimeSpentInA_before = cell(numel(genAll),1);
percentTimeSpentInA_during = cell(numel(genAll),1);

ttInA = @(x) sum(x<=1);
ttInB = @(x) sum(x>1 & x<=1.3);
for gen = 1:numel(genAll)
    
    load([meta.folderObject  '\' genAll{gen} '_' meta.d meta.ext '.mat'],'f_orco')
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

figure;set(gcf,'Position',[2 42 838 924]);
k = 1;
opts.plotData = false;
opts.yl = [0 1];
for gen = 1:numel(genAll)
    for gen2 = (gen+1):numel(genAll)
        data = nan(max([nFly(gen) nFly(gen2)]),2);
        data(1:nFly(gen),1) = percentTimeSpentInA_during{gen};
        data(1:nFly(gen2),2) = percentTimeSpentInA_during{gen2};
        opts.xLabels = genAll([gen gen2]);
        subplot(3,2,k)
        violinPlotsStats(data,opts)
        xtickangle(10)
        k = k+1;
    end
end

