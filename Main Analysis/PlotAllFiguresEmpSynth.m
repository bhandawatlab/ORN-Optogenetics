function [] = PlotAllFiguresEmpSynth(genAll,meta,agentModelMeta)

% set parameters
tau = agentModelMeta.tau;
C = agentModelMeta.C;

plotFolder = [meta.syntheticPlotFold];

for ii = 1:numel(tau)
    for jj = 1:numel(C)
        currTau = tau(ii);
        currC = C(jj);
        cond = ['_CWTau' strrep(num2str(currTau),'.','') ...
            '_CWC' strrep(num2str(currC),'.','')];
        
        f_orcoAll = cell(1,numel(genAll));
        synth_orcoAll = cell(1,numel(genAll));
        k = 1;tic;
        for g = 1:numel(genAll)
            gen = genAll{g};
            fName = ['RT_run' gen '_' meta.d meta.ext cond '_flies.mat'];
            load([meta.syntheticFlyFold fName],'f_orco','synth_orco')
            
            f_orco.id = [f_orco.id ' tau=' num2str(currTau) ', C=' num2str(currC)];
            synth_orco.id = [synth_orco.id ' tau=' num2str(currTau) ', C=' num2str(currC)];
            f_orcoAll{k} = f_orco;
            synth_orcoAll{k} = synth_orco;
            k = k+1;
        end
        fprintf('Finished loading in data in %d seconds\n',round(toc))
        
        ext = '';%['EnterLeaveExpDecay_' gen];
        
        %%
        % plot radial occupancy
        close all;fNum = 1;tic;
        stopSpd = meta.stopThresh;
        [fNum, ~] = plotRadialOccupancy(f_orcoAll,synth_orcoAll,meta.border,stopSpd,fNum);
        printFigures(fNum,plotFolder,['EmpSynth_RadialOccupancy_noStops' ext])
        fprintf('Finished plotting radial occupancy in %d seconds\n',round(toc))
        
        %%
        % plot radial occupancy
        close all;fNum = 1;tic;
        stopSpd = 0;%meta.stopThresh;
        [fNum, radProbAllGen] = plotRadialOccupancy(f_orcoAll,synth_orcoAll,meta.border,stopSpd,fNum);
        printFigures(fNum,plotFolder,['EmpSynth_RadialOccupancy' ext])
        fprintf('Finished plotting radial occupancy in %d seconds\n',round(toc))

        %%
        % plot turn density
        close all;fNum = 1;tic;
        [fNum, turnDensAllGen] = plotTurnDensity(f_orcoAll,synth_orcoAll,meta.border,fNum);
        %save(['allGenModelRadPosTurnDens' meta.ext '.mat'],'radProbAllGen','turnDensAllGen','genAll','-v7.3');
        printFigures(fNum,plotFolder,['EmpSynth_turnDensity2' ext])
        fprintf('Finished plotting turn density in %d seconds\n',round(toc))
        
        %%
        % plot correlation between turn density/radial occupancy and model
        % fit
        close all;fNum = 1;tic;
        [fNum] = plotTurnDensityRadialOccupencyCorr(radProbAllGen,turnDensAllGen,genAll,fNum);
        printFigures(fNum,plotFolder,['EmpSynth_correlationAnalysis' ext])
        fprintf('Finished plotting correlation analysis in %d seconds\n',round(toc))
        
        %%
        % plot probability of inside
        close all;fNum = 1;tic;
        fNum = plotProbInside(f_orcoAll,synth_orcoAll,meta.border,fNum);
        printFigures(fNum,plotFolder,['EmpSynth_ProbInside' ext])
        fprintf('Finished plotting probability of being inside in %d seconds\n',round(toc))
        
        %%
        % plot attraction index (before vs during)
        close all;fNum = 1;tic;
        fNum = plotAttnNdxBefDur(f_orcoAll,synth_orcoAll,meta.border,fNum);
        % plot attraction index (emp vs synth)
        fNum = plotAttnNdxEmpSyth(f_orcoAll,synth_orcoAll,meta.border,fNum);
        printFigures(fNum,plotFolder,['EmpSynth_AttractionNdx' ext])
        fprintf('Finished plotting attraction index in %d seconds\n',round(toc))
        
        %%
%         % plot spatialtemporal density
%         close all;fNum = 1;tic
%         fNum = plotSpatialTemporalDensity(f_orcoAll,synth_orcoAll,fNum);
%         printFigures(fNum,plotFolder,['EmpSynth_SpatialTemporal_Density' ext])
%         fprintf('Finished spatial-temporal density in %d seconds\n',round(toc))
        
        %%
        % plot trajectories
        close all;tic;
        for g = 1:numel(genAll)
            close all;fNum = 1;
            f_orco = f_orcoAll{g};
            synth_orco = synth_orcoAll{g};
            plotXYTracks(f_orco,meta.border,[],true,'Empirical')
            plotXYTracks(synth_orco,meta.border,[],true,'Synthetic')
            fNum = get(gcf,'Number')+1;
            printFigures(fNum,[plotFolder 'XY tracks/'],['EmpSynth_' f_orco.id '_tracksByAttaction'])
        end
        fprintf('Finished plotting trajectories in %d seconds\n',round(toc))
        
        
    end
end
end

function [fNum] = plotAttnNdxBefDur(f_orcoAll,synth_orcoAll,border,fNum)
k = 1;cond = {'Empirical','Synthetic'};
for g = 1:numel(f_orcoAll)
    currEmpSynth = {f_orcoAll{g},synth_orcoAll{g}};
    
    if mod(k,24) == 1
        subplots = [6,4,1,fNum]; % ysbplt, xsbplt, sbpltN, fig
        xSubplots = repmat(100./subplots(1),1,subplots(1));
        ySubplots = repmat(100./subplots(2),1,subplots(2));
        
        figure(subplots(4));set(gcf,'Position',[2 42 838 924])
        p = panel();
        p.pack(xSubplots, ySubplots);
        subplots(3) = 1;
        
        % define basic paramters
        lims = [0 1];                   % y limits to set
        isPaired = 'N';                 % equivalent to paired t-test
        circleSize = 60;                % size of scatter plot points
        barstate = 'off';               % 'On' = bar graph, 'Off' = scatter plot
        
        fNum = fNum+1;
    end
    
    for i = 1:2
        attNdx = currEmpSynth{i}.getAttnNdx('H',border);
        
        dat1 = [attNdx.before attNdx.during];
        g1 = [repmat({'before'},currEmpSynth{i}.nFly,1); repmat({'during'},currEmpSynth{i}.nFly,1)];
        [ss,p,~] = dabest3(dat1,g1,p,[],lims,isPaired,circleSize,barstate,subplots);
        % set axis labels for comparisons against category 1
        [J,I] = ind2sub(subplots(2:-1:1),subplots(3));
        try
            p(I,J).select();
            text(0.2,0.8,['n= ' num2str(numel(dat1)./2) ' flies'],'Color','k')
            text(0.2,0.9,cond{i},'Color','k')
            ylabel('atn Ndx');
            title(currEmpSynth{i}.id,'Interpreter','none')
        catch
            p(I,J,1,1).select();
            ylabel('atn Ndx');
            title(currEmpSynth{i}.id,'Interpreter','none')
            p(I,J,2,1).select();
            ylabel('delta atn Ndx');
            xtickangle(10)
        end
        subplots(3) = subplots(3)+1;
    end
    k = subplots(3);
end
end

function [fNum] = plotAttnNdxEmpSyth(f_orcoAll,synth_orcoAll,border,fNum)

k = 1;cond = {'Emp','Synth'};scenario = {'Before','During'};
for g = 1:numel(f_orcoAll)
    currEmpSynth = {f_orcoAll{g},synth_orcoAll{g}};
    
    if mod(k,24) == 1
        subplots = [6,4,1,fNum]; % ysbplt, xsbplt, sbpltN, fig
        xSubplots = repmat(100./subplots(1),1,subplots(1));
        ySubplots = repmat(100./subplots(2),1,subplots(2));
        
        figure(subplots(4));set(gcf,'Position',[2 42 838 924])
        p = panel();
        p.pack(xSubplots, ySubplots);
        subplots(3) = 1;
        
        % define basic paramters
        lims = [0 1];                   % y limits to set
        isPaired = 'N';                 % equivalent to paired t-test
        circleSize = 60;                % size of scatter plot points
        barstate = 'off';               % 'On' = bar graph, 'Off' = scatter plot
        
        fNum = fNum+1;
    end
    
    attNdx_befDur = cell(1,2);
    g1 = {};n = zeros(1,2);
    for i = 1:2
        attNdx = currEmpSynth{i}.getAttnNdx('H',border);
        
        attNdx_befDur{1} = [attNdx_befDur{1} attNdx.before];
        attNdx_befDur{2} = [attNdx_befDur{2} attNdx.during];
        
        g1 = [g1;repmat(cond(i),currEmpSynth{i}.nFly,1)];
        n(i) = numel(attNdx.before);
    end
    
    for i = 1:2
        dat1 = attNdx_befDur{i};
        [ss,p,~] = dabest3(dat1,g1,p,[],lims,isPaired,circleSize,barstate,subplots);
        % set axis labels for comparisons against category 1
        [J,I] = ind2sub(subplots(2:-1:1),subplots(3));
        try
            p(I,J).select();
            text(0.2,0.7,['n= ' num2str(n(2)) ' synth flies'],'Color','k')
            text(0.2,0.8,['n= ' num2str(n(1)) ' emp flies'],'Color','k')
            text(0.2,0.9,scenario{i},'Color','k')
            ylabel('atn Ndx');
            title(currEmpSynth{i}.id,'Interpreter','none')
        catch
            p(I,J,1,1).select();
            ylabel('atn Ndx');
            title(currEmpSynth{i}.id,'Interpreter','none')
            p(I,J,2,1).select();
            ylabel('delta atn Ndx');
            xtickangle(10)
        end
        subplots(3) = subplots(3)+1;
    end
    k = subplots(3);
end
end

function [fNum] = plotProbInside(f_orcoAll,synth_orcoAll,border,fNum)

cc = {'k','r'};k = 1;
for g = 1:numel(f_orcoAll)
    currEmpSynth = {f_orcoAll{g},synth_orcoAll{g}};
    
    if mod(k,24) == 1
        figure(fNum);set(gcf,'Position',[2 42 838 924])
        suptitle('Probability Inside')
        k = 1;fNum = fNum+1;
    end
    
    for i = 1:2
        fe = currEmpSynth{i}.getFirstEntry('H',border);
        probIn = currEmpSynth{i}.getProbIn('H',border,fe);
        probIn2 = currEmpSynth{i}.getProbIn('H',border,ones(currEmpSynth{i}.nFly,1));
        
        % prob being inside %
        subplot(6,4,k);
        shadedErrorBar([1:1:numel(probIn2.mean)]./currEmpSynth{i}.fs./60,...
            probIn2.mean,probIn2.sem,'lineprops',cc{i});hold on;% in minutes
        
        subplot(6,4,k+1);
        shadedErrorBar([1:1:numel(probIn.mean)]./currEmpSynth{i}.fs./60,...
            probIn.mean,probIn.sem,'lineprops',cc{i});hold on;% in minutes
    end
    
    subplot(6,4,k);
    xlim([0 6]);ylim([0 0.6])
    xlabel('time (min)')
    ylabel('Probability')
    title(currEmpSynth{i}.id)
    
    subplot(6,4,k+1);
    xlim([0 3]);ylim([0 1])
    xlabel('time since first entry (min)')
    ylabel('Probability')
    title(currEmpSynth{i}.id)
    
    if k==1
        subplot(6,4,k);title({currEmpSynth{i}.id, ' black=emp,red=synth'})
    end
    
    k = k+2;
end
end

function [fNum,radProbAllGen] = plotRadialOccupancy(f_orcoAll,synth_orcoAll,border,stopSpd,fNum)
cc = {'k','r'};k = 1;
radProbAllGen = cell(1,numel(f_orcoAll));
for g = 1:numel(f_orcoAll)
    currEmpSynth = {f_orcoAll{g},synth_orcoAll{g}};
    
    if mod(k,24) == 1
        figure(fNum);set(gcf,'Position',[2 42 838 924])
        suptitle('Radial Occupancy')
        k = 1;fNum = fNum+1;
    end
    dur = 0;
    for i = 1:2
        %radProb = currEmpSynth{i}.getRadialProb('H',border,0,currEmpSynth{i}.fs);
        radProb = currEmpSynth{i}.getRadialProb('H',border,stopSpd,dur);
        
        subplot(6,4,k);
        % plot radial probability
        plot(radProb.before.x.*currEmpSynth{i}.rBound,...
            radProb.before.y,'Color',cc{i},'Linewidth',1);hold on;
        
        subplot(6,4,k+1);
        plot(radProb.during.x.*currEmpSynth{i}.rBound,...
            radProb.during.y,'Color',cc{i},'Linewidth',1);hold on;
        
        
        radProbAllGen{i,g} = radProb;
    end
    
    subplot(6,4,k);
    plot([border,border],[0 0.4],'r--');ylim([0 0.3])
    xlabel('Position (cm)');ylabel('Probability');
    title({currEmpSynth{i}.id, ' Before'})
    
    subplot(6,4,k+1);
    plot([border,border],[0 0.4],'r--');ylim([0 0.3])
    xlabel('Position (cm)');ylabel('Probability');
    title({currEmpSynth{i}.id, ' During'})
    
    if k==1
        subplot(6,4,k);legend({'emp Head','synth'})
    end
    
    k = k+2;
end
end

function [fNum,turnDensAllGen] = plotTurnDensity(f_orcoAll,synth_orcoAll,border,fNum)

cc = {'k','r'};k = 1;
turnDensAllGen = cell(1,numel(f_orcoAll));
for g = 1:numel(f_orcoAll)
    f_orco = f_orcoAll{g};
    synth_orco = synth_orcoAll{g};
    
    currEmpSynth = {f_orco,synth_orco};
    
    if mod(k,24) == 1
        figure(fNum);set(gcf,'Position',[2 42 838 924])
        suptitle('Radial turn density')
        k = 1;fNum = fNum+1;
    end
    
    for type = 1:2
        % radial location of turns
        %-------------
        fe = currEmpSynth{type}.getFirstEntry('H',border);
        spk = currEmpSynth{type}.spk;
        for j = 1:currEmpSynth{type}.nFly
            spk(j,:) = smooth(spk(j,:));
        end
        dSpk = gradient(spk);
        baseline = f_orco.spk(1);
        
        condNdx = zeros(size(spk));condNdx(spk>baseline) = 4;
        condNdx(spk<baseline) = 2;condNdx(abs((spk-baseline))<0.001 & dSpk==0) = 3;
        for j = 1:currEmpSynth{type}.nFly
            condNdx(j,1:fe(j)-1) = 1;
        end
        thresh = 0.1;
        key = {'before','below','baseline','above'};
        [s,keyNdx] = currEmpSynth{type}.getKinematicsCond(condNdx,key,thresh,[0 200],true);%%
        turnDens = currEmpSynth{type}.getTurnDens('H',s,keyNdx);
        
        turnDens.x(end) = 1;
        xx = (turnDens.x(1:end-1)+turnDens.x(2:end))./2;
        
        subplot(6,4,k);
        plot(xx,turnDens.Bef,cc{type},'Linewidth',2);hold on
        
        subplot(6,4,k+1);
        plot(xx,turnDens.Dur,cc{type},'Linewidth',2);hold on;
        
        turnDensAllGen{type,g} = turnDens;
    end
    
    subplot(6,4,k);
    ylim([0 0.2])
    xticks([0:0.5:4]./4)
    title({currEmpSynth{1}.id, ' Before'})
    
    subplot(6,4,k+1);
    ylim([0 0.2])
    xticks([0:0.5:4]./4)
    title({currEmpSynth{1}.id, ' During'})
    
    if k == 1
        %legend({'emp Head','synth'})
    end
    
    k = k+2;
end

end

function [fNum] = plotSpatialTemporalDensity(f_orcoAll,synth_orcoAll,fNum)

k = 1;
for g = 1:numel(f_orcoAll)
    fprintf('Generating RadialPositionPlots, %d / %d\n',[g, numel(f_orcoAll)]);
    f_orco = f_orcoAll{g};
    synth_orco = synth_orcoAll{g};
    [y_emp,xCentEmp,yCentEmp] = plotRadialPosition(f_orco,'after',false,'dt',2.*f_orco.fs,'clims',[0 0.01]);
    %[y_synth,xCent,yCent] = plotRadialPosition(f_orco,'after',false,'dt',2.*f_orco.fs);toc;
    [y_synth,xCentSynth,yCentSynth] = plotRadialPosition(synth_orco,'after',false,'dt',2.*synth_orco.fs);
    
    if mod(k,24) == 1
        figure(fNum);set(gcf,'Position',[2 42 838 924])
        k = 1;fNum = fNum+1;
    end
    clims = [0 ceil(max([y_emp(:);y_synth(:)]).*1000)./1000];
    subplot(6,4,k);
    imagesc(xCentEmp./f_orco.fs, yCentEmp./4,y_emp,clims);colorbar
    xlabel('time');ylabel('r')
    shading interp
    title({f_orco.id 'Emp'})
    xlim([0 180]);ylim([0 4]);
    xticks([0:60:180]);xticklabels(cellstr(num2str([0:60:180]')))
    yticks([0:2:4]); yticklabels(cellstr(num2str([0:2:4]')))
    
    subplot(6,4,k+1);
    imagesc(xCentSynth./synth_orco.fs, yCentSynth./4,y_synth,clims);colorbar
    xlabel('time');ylabel('r')
    shading interp
    title({f_orco.id 'Synth'})
    xlim([0 180]);ylim([0 4]);
    xticks([0:60:180]);xticklabels(cellstr(num2str([0:60:180]')))
    yticks([0:2:4]); yticklabels(cellstr(num2str([0:2:4]')))
    
    k = k+2;
end
end

function [] = printFigures(fNum,folder,figTitle)
if ~exist(folder, 'dir')
    mkdir(folder)
end
if numel(fNum)==1
    fRange = 1:fNum-1;
else
    fRange = fNum;
end
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
