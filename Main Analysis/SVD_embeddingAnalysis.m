function [] = SVD_embeddingAnalysis(f_orco,meta,figureFile)
% SVD_embeddingAnalysis  Compares state space dimensions from taking 
%   derivatives vs from truncated SVD. Uses the false nearest neighbor
%   algorithm to determine embedding dimensions
%
%   Inputs: f_orco = a fly object
%           meta = a structure with field border for light border radius
%           and field plotFig for whether to save figure
%           figureFile = output figure name is meta.plotFig = true
%   

%load('DataModel\Orco Retinal_April2022_allTime.mat', 'f_orco')
fe = f_orco.getFirstEntry('H',meta.border);
opts.ncells = 1;

%%
close all
eLag = 50;
FNN = 100.*ones(f_orco.nFly,5);FNN_svd = 100.*ones(f_orco.nFly,5);
for fly = 1:f_orco.nFly
    % get the 0th - 4th derivatives of the firing rate
    xdata = f_orco.spk(fly,(fe(fly)-30):end)';
    stateData = smooth(xdata+(rand(size(xdata))-0.5).*0.001,6);
    for d = 2:5
        stateData(:,d) = gradient(stateData(:,d-1));
    end
    
    % calculate the SVD of the hankel matrix
    ntfilt = 6;% choose a 200 ms delay embedding
    [XStimAll,~,~,~] = createDesignMat(xdata',[],ntfilt,[],opts);
    XStimAll = XStimAll(eLag:end,2:end);
    [U,~,~] = svd(XStimAll);
    
    % calculate False nearest neighbor for dims based on derivative
    for d = 2:5
        xdata_d = stateData(:,1:d-1);
        xdata_d_p = stateData(:,1:d);
        [Idx,D] = knnsearch(xdata_d,xdata_d,'K',2,'Distance','seuclidean');
        D = D(:,2);
        D2 = sum((xdata_d_p(Idx(:,2),:)-xdata_d_p).^2,2);
        FNN(fly,d) = sum(sqrt((D2-D)./D)>10)./size(D,1).*100;
        %D2 = D+D2(:,2:end);
    end
    
    % calculate False nearest neighbor for dims based on truncated SVD
    for d = 2:5
        xdata_d = U(:,1:d-1);
        xdata_d_p = U(:,1:d);
        [Idx,D] = knnsearch(xdata_d,xdata_d,'K',2,'Distance','seuclidean');
        D = D(:,2);
        D2 = sum((xdata_d_p(Idx(:,2),:)-xdata_d_p).^2,2);
        FNN_svd(fly,d) = sum(sqrt((D2-D)./D)>10)./size(D,1).*100;
    end
    
    % set up figure
    if mod(fly,3)==1
        figure;set(gcf,'Position',[2 42 838 924]);
        k = 1;
    end
    
    % plot the 0,1,2nd derivative of the firing rate and the first 3
    % dimensions from tSVD
    subplot(3,2,k);
    plot3(-stateData(eLag:end,1),stateData(eLag:end,2),stateData(eLag:end,3));view([225 45]);%view([45 45])
    xlabel('-F');ylabel('dF');zlabel('d(dF)')
    title(['Fly' num2str(fly) ' raw'])
    subplot(3,2,k+1);
    plot3(U(:,1),U(:,2),U(:,3));view([225 45]);%view([45 45])
    xlabel('dim 1');ylabel('dim 2');zlabel('dim 3')
    title(['Fly' num2str(fly) ' SVD'])
    k = k+2;
end

% plot the false nearest neighbor results
figure;set(gcf,'Position',[2 42 838 924]);
subplot(2,1,1);
shadedErrorBar(1:5,nanmean(FNN),std(FNN));hold on;
h = plot((find((nanmean(FNN)<1),1)+1).*[1 1],[0 100],'--r');hold off;
ylim([0 100])
legend(h,'embedding dimensions')
xlabel('dimensions');ylabel('Percent FNN');
title('derivatives')
subplot(2,1,2);
shadedErrorBar(1:5,nanmean(FNN_svd),std(FNN_svd));hold on;
h = plot((find((nanmean(FNN_svd)<1),1)+1).*[1 1],[0 100],'--r');hold off;
ylim([0 100])
legend(h,'embedding dimensions')
xlabel('dimensions');ylabel('Percent FNN');
title('SVD')

if meta.plotFig
    for f = 1:get(gcf,'Number')
        figure(f);
        print('-painters','-dpsc2',[figureFile '.ps'],'-loose','-append');
    end
    ps2pdf('psfile', [figureFile '.ps'], 'pdffile', ...
    [figureFile '.pdf'], 'gspapersize', 'letter',...
    'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
    'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
    'gslibpath','C:\Program Files\gs\gs9.50\lib');
end


end