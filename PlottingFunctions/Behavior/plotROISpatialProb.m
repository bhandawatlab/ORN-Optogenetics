function [fNum] = plotROISpatialProb(f_orcoAll,dxy,lab,border,fNum)
if isempty(lab)
    lab = '';
end

T_df = [-20,20];
T_f = 15;

ROI_func{6} = @(f,df,baseline) (f-baseline)<0.1 & (df)<0.1;% baseline
ROI_func{1} = @(f,df,baseline) df<T_df(1) & ~ROI_func{6}(f,df,baseline);% negative df
ROI_func{2} = @(f,df,baseline) df>T_df(2) & ~ROI_func{6}(f,df,baseline);% positive df
ROI_func{3} = @(f,df,baseline) f>T_f & ~ROI_func{1}(f,df,baseline) & ...
    ~ROI_func{2}(f,df,baseline) & ~ROI_func{6}(f,df,baseline);% high f
ROI_func{4} = @(f,df,baseline) f<baseline & ~ROI_func{1}(f,df,baseline) & ...
    ~ROI_func{2}(f,df,baseline) & ~ROI_func{6}(f,df,baseline);% below baseline
ROI_func{5} = @(f,df,baseline) ~ROI_func{1}(f,df,baseline) & ...
    ~ROI_func{2}(f,df,baseline) & ~ROI_func{3}(f,df,baseline) & ...
    ~ROI_func{4}(f,df,baseline) & ~ROI_func{6}(f,df,baseline);% low f
ROI_label = {'Neg delta F','Pos delta F','High F','Inhibition','Low F','baseline'};

k = 1;
nGen = numel(f_orcoAll);
for gen = 1:nGen
    f_orco = f_orcoAll{gen};
    fe = f_orco.getFirstEntry('H',border);
    
    fe_mask = false(size(f_orco.x));
    for fly = 1:numel(fe)
        fe_mask(fly,fe(fly):end) = true;
    end
    
    hist = ceil(0.2.*f_orco.fs)-1;%200 ms history
    
    dSpkRateAll = f_orco.calcDeltaFR;
    
    dFSmooth = dSpkRateAll(:,1:end-hist);
    FSmooth = f_orco.spk(:,1:end-hist);
    for i = 2:hist+1
        dFSmooth = dSpkRateAll(:,i:end-(hist+1)+i)+dFSmooth;
        FSmooth = f_orco.spk(:,i:end-(hist+1)+i)+FSmooth;
    end
    dFSmooth = [zeros(size(dFSmooth,1),hist),dFSmooth./(hist+1)];% average
    FSmooth = [f_orco.spk(1).*ones(size(FSmooth,1),hist),FSmooth./(hist+1)];% average
    
    baseline = FSmooth(1);
        
    xPos = f_orco.(['x' lab]);
    yPos = f_orco.(['y' lab]);
    rPos = f_orco.(['r' lab]);
    
    xPos = xPos(fe_mask);
    yPos = yPos(fe_mask);
    rPos = rPos(fe_mask);
    f = FSmooth(fe_mask);
    df = dFSmooth(fe_mask);
    
    % section by x,y position
    Xedges = -f_orco.rBound:dxy(1):f_orco.rBound;
    Yedges = -f_orco.rBound:dxy(2):f_orco.rBound;
    [N,Xedges,Yedges,binX,binY] = histcounts2(xPos,yPos,Xedges,Yedges);
    
    % section by r position
    Redges = 0:dxy(1):f_orco.rBound;
    [N_r,Redges,binR] = histcounts(rPos,Redges);
    
    % get mask for each region
    for i = 1:6
        ROI_mask{i} = ROI_func{i}(f,df,baseline);
    end
    
    % count the number of data points for each region given x,y position
    tic;
    N_ROI = zeros([size(N),6]);
    for x = 1:size(N,1)
        for y = 1:size(N,2)
            c = binX==x & binY==y;
            for i = 1:6
                N_ROI(x,y,i) = sum(c & ROI_mask{i},'all');
            end
        end
    end
    toc;
    
    % count the number of data points for each region given radial position
    tic;
    N_r_ROI = zeros([size(N_r,2),6]);
    for r = 1:size(N_r,2)
        c = binR==r;
        for i = 1:6
            N_r_ROI(r,i) = sum(c & ROI_mask{i},'all');
        end
    end
    toc;
    
    if gen == 1
        N_ROI_All = N_ROI;
        N_All = N;
        N_r_ROI_All = N_r_ROI;
        N_r_All = N_r;
    else
        N_ROI_All = N_ROI_All+N_ROI;
        N_All = N_All+N;
        N_r_ROI_All = N_r_ROI_All+N_r_ROI;
        N_r_All = N_r_All+N_r;
    end
    
    xBin = (Xedges(1:end-1)+Xedges(2:end))./2;
    yBin = (Yedges(1:end-1)+Yedges(2:end))./2;
    rBin = (Redges(1:end-1)+Redges(2:end))./2;
    
    for i = 1:6
        if mod(k,96) == 1
            figure(fNum);set(gcf,'Position',[2 42 838 924])
            suptitle({'spatial grid: # ROI/# total data',' ',' '})
            k = 1;fNum = fNum+1;
        end
        
        subplot(16,6,[k k+6]);
        imagesc(xBin,yBin,N_ROI(:,:,i)./N,[0 1]);
        colormap([0,0,0;jet]);
        if i == 1
            title({f_orco.id; ROI_label{i}})
        else
            title({' '; ROI_label{i}})
        end
        k = k+1;
    end
    k = k+6;

    if mod(k+11,96)==0
        plotCB([0:0.2:1])
    end
    
    for i = 1:6
        if mod(k,96) == 1
            figure(fNum);set(gcf,'Position',[2 42 838 924])
            suptitle({'spatial grid: # ROI/# total data',' ',' '})
            k = 1;fNum = fNum+1;
        end
        tmp = [flipud(N_r_ROI(:,i)); N_r_ROI(:,i)]'./[fliplr(N_r),N_r];
        subplot(16,6,k);
        plot([-fliplr(rBin),rBin],tmp)
        xlim([-f_orco.rBound f_orco.rBound]);ylim([0 ceil(max(tmp).*10)./10])
        k = k+1;
    end
    k = k+6;
end

for i = 1:6
    if mod(k,96) == 1
        figure(fNum);set(gcf,'Position',[2 42 838 924])
        suptitle({'spatial grid: # ROI/# total data',' ',' '})
        k = 1;fNum = fNum+1;
    end
    
    subplot(16,6,[k k+6]);
    imagesc(xBin,yBin,N_ROI_All(:,:,i)./N_All,[0 1]);
    colormap([0,0,0;jet]);
    title({'All genotypes'; ROI_label{i}})
    k = k+1;
end
plotCB([0:0.2:1])
k = k+6;

for i = 1:6
    if mod(k,96) == 1
        figure(fNum);set(gcf,'Position',[2 42 838 924])
        suptitle({'spatial grid: # ROI/# total data',' ',' '})
        k = 1;fNum = fNum+1;
    end
    tmp = [flipud(N_r_ROI_All(:,i)); N_r_ROI_All(:,i)]'./[fliplr(N_r_All),N_r_All];
    subplot(16,6,k);
    plot([-fliplr(rBin),rBin],tmp)
    xlim([-f_orco.rBound f_orco.rBound]);ylim([0 ceil(max(tmp).*10)./10])
    k = k+1;
end

end


function [] = plotCB(ztick)

cb=colorbar('XTick', ztick);
cb.Position = cb.Position + [0.45e-1, 0.1e-2, 0, 0];

end