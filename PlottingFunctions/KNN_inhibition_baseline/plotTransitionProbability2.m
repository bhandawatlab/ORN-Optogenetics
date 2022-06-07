function [fNum] = plotTransitionProbability2(f_orcoAll,plotType,colortype,fNum)

% loop through each genotype
k=1;
nGen = numel(f_orcoAll);
for gen = 1:nGen
    f_orco = f_orcoAll{gen};
    m = f_orco.model.TP.during;
    XX = f_orco.model.TP.XX;
    YY = f_orco.model.TP.YY;
    
    nStatePrior = numel(f_orco.model.TP.key);
    nStateNew = nStatePrior-1;
    
    for i = 1:nStatePrior
        priorState = f_orco.model.TP.key{i};
        for j = 1:nStateNew
            newState = f_orco.model.TP.key{j};
            m = f_orco.model.TP.during{i,j};
            mBaseline = f_orco.model.TP.before(i,j);
            if ~isempty(m)
                ylims = [0 1];ztick = [0:0.5:1];

                % average across time
                if mod(k,24) == 1
                    figure(fNum);set(gcf,'Position',[2 42 838 924])
                    currSupT = [f_orcoAll{gen}.id];
                    for g = 1:min(nGen-gen,2)
                        currSupT = [currSupT '; ' f_orcoAll{gen+g}.id];
                    end
                    currSupT = [currSupT '; time from first entry'];
                    suptitle(currSupT)
                    k = 1;fNum = fNum+1;
                end

                subplot(6,4,k);
                plotSubplot(XX, YY, m, mBaseline, ylims,ztick, colortype, plotType);
                title([priorState ' to ' newState])
                plotCB(mBaseline,ztick)
                if mod(k,4)~=1
                    ylabel('')
                end
                k = k+1;
            end
        end
    end
    k = ceil((k-1)/4)*4+1;
    if k>24
        k = 25;
    end
end


end

function [] = plotSubplot(XX, YY, m, mBaseline, ylims, ztick, colortype, plotType)
switch colortype
    case 'absolute'
        ylims(1) = ylims(1);
    case 'relative'
        baseline = mBaseline;
        if baseline<=mean(ylims)
            ylims(1) = baseline-(ylims(2)-baseline);
        else
            ylims(2) = baseline.*2;
        end
end

switch plotType
    case 'imagesc'
        imagesc(XX(1,:),YY(:,1)',m(XX,YY),[ylims(1) ylims(2)]);set(gca,'YDir','normal');
        colormap([0 0 0;jet]);
        %plotCB(m,ztick)
    case 'surf'
        surf(XX,YY,exp(m.val(:,:,tt)));view([45 45]);hold on;
        points=[[min(XX(:)) min(XX(:)) max(XX(:)) max(XX(:))]' ...
            [0 max(YY(:)) max(YY(:)) 0]' mBaseline.*[1 1 1 1]'];
        h = fill3(points(:,1),points(:,2),points(:,3),'k');
        grid on
        alpha(h, 0.1)
        
        caxis([ylims(1) ylims(2)]);
        colormap(jet)
        zlim(ylims);zticks(ztick);
        colormap([0 0 0;jet]);
        hold off;
        %plotCB(m,ztick)
end

xlabel('dF (spks/s^2)');ylabel('F (spks/s)');
xlim([-150 150]);ylim([0 50]);
xticks([-150:75:150]);yticks([0:25:50]);

end

function [] = plotCB(mBaseline,ztick)

cb=colorbar('XTick', ztick);cb.Position = cb.Position + [0.45e-1, 0.1e-2, 0, 0];
[newTicks,ndx] = sort([cb.Ticks, mBaseline-eps]);
cb.TickLabels = cellstr(num2str(cb.Ticks'));
newLabels = [cb.TickLabels; 'baseline'];
cb.Ticks = newTicks;
cb.TickLabels = newLabels(ndx);

end