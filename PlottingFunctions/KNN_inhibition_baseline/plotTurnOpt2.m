function [fNum,stateFigure] = plotTurnOpt2(f_orcoAll,state2Cons,tim2Cons,plotType,...
    colortype,plotAdaptation,plotInterp,fNum)

k=ones(1,numel(state2Cons));
nGen_per_page = 3;
nPage_all_gen = ceil(numel(f_orcoAll)./nGen_per_page);
fNum_byState = fNum:nPage_all_gen:numel(state2Cons).*nPage_all_gen;
nGen = numel(f_orcoAll);
for gen = 1:nGen
    f_orco = f_orcoAll{gen};
    %m = f_orco.model.TurnBias;
    
    % turn optimality
    m = f_orco.model.TurnBias;
    tGrid = m.tt;
    ztick = [0:0.25:1];
    for state = 1:numel(state2Cons)
        
        % average across time
        if numel(tim2Cons)==1
            if mod(k(1),12) == 1
                figure(fNum);set(gcf,'Position',[2 42 838 924])
                currSupT = [f_orcoAll{gen}.id];
                for g = 1:min(nGen-gen,3)
                    currSupT = [currSupT '; ' f_orcoAll{gen+g}.id];
                end
                currSupT = {currSupT, ' Turn Optimality time from first entry'};
                sgtitle(currSupT)
                
                stateFigure{fNum} = 'all';
                k(1) = 1;fNum = fNum+1;
            end
            tt = tim2Cons(1)==tGrid;
            
            subplot(4,3,k(1));
            ylims = plotSubplot(m, state, tt, tGrid, colortype, plotType, plotInterp);
            plotCB(m,ztick,ylims,state)
            if mod(k(1),3)~=1
                ylabel('')
            end
            k(1) = k(1)+1;
            
            [XAll{gen,state,tt},yAll{gen,state,tt},OLS{gen,state}(:,tt),...
                OLSdFOnly{gen,state}(:,tt)] = calcOLS(m,tt,state);
        else
            if mod(k(state),24) == 1
                figure(fNum_byState(state));set(gcf,'Position',[2 42 838 924])
                currSupT = [f_orcoAll{gen}.id];
                for g = 1:min(nGen-gen,2)
                    currSupT = [currSupT '; ' f_orcoAll{gen+g}.id];
                end
                currSupT = [currSupT ' Turn Optimality time from first entry'];
                suptitle(currSupT)
                stateFigure{fNum_byState(state)} = m.key{state};
                
                k(state) = 1;fNum_byState(state) = fNum_byState(state)+1;
            end
            figure(fNum_byState(state)-1);
            
            for t = 1:min(numel(tim2Cons),8)%24
                tt = tim2Cons(t)==tGrid;
                subplot(6,4,k(state));
                ylims = plotSubplot(m, state, tt, tGrid, colortype, plotType, plotInterp);
                k(state)=k(state)+1;
            end
            for tt = 1:min(numel(tGrid),24)
                [XAll{gen,state,tt},yAll{gen,state,tt},OLS{gen,state}(:,tt),...
                    OLSdFOnly{gen,state}(:,tt)] = calcOLS(m,tt,state);
            end
            k(state) = k(state)+(8-t);
            
            if mod(k(state)-1,8)==0
                plotCB(m,ztick,ylims,state)
            end
        end
    end
    
    if numel(tim2Cons)==1
        k(1) = k(1)+(3-numel(state2Cons));
    end
    
end

fNum = max(fNum,fNum_byState(end));

if plotAdaptation
    figure(fNum);set(gcf,'Position',[2 42 838 924])
    for state = 1:3
        subplot(3,1,state)
        scatter(XAll{1,state,1}(:,2),yAll{1,state,1});hold on;
        plot([-150:150],[ones(1,301);-150:150]'*OLSdFOnly{1,state}(:,1))
        ylim([0 1]);xlim([-150 150])
        ylabel('Probability');xlabel('dF (spks/s^2)');
        legend({'all KNN grids','OLS fit'})
        title([m.key{state} ' optimality'])
    end
    suptitle('Example OLS fit')
    fNum = fNum+1;
    
    tt = tGrid(1:min(numel(tGrid),24));
    k = 1;
    for gen = 1:numel(f_orcoAll)
        f_orco = f_orcoAll{gen};
        if mod(gen,12) == 1
            figure(fNum);set(gcf,'Position',[2 42 838 924])
            suptitle('Turn optimality dF slope time from first entry')
            fNum = fNum+1;
            k = 1;
        end
        subplot(4,3,k);hold on
        for state = 1:3
            plot(tt,OLSdFOnly{gen,state}(2,:),'Linewidth',2);
        end
        plot([0 tt(end)],[0 0],'--k','Linewidth',2);hold off;
        if mod(k,12)==1
            legend([f_orco.states.key(1:3)])
        end
        xlim([0 max(tt)]);ylim([-1 1]./300)
        xticks([0:30:max(tt)]);yticks((-3:3).*10.^-3);
        title([f_orco.id])
        xlabel('time since first entry (s)');ylabel('dOpt/dF (1/(spk/s^2))')
        k=k+1;
    end
end

end

function [ylims] = plotSubplot( m, state, tt, tGrid, colortype, plotType, plotInterp)
switch colortype
    case 'absolute'
        ylims = [0 1];
    case 'relative'
        ylims = [0 1];
        ylims(1) = m.before{state}-(ylims(2)-m.before{state});
end
if plotInterp
    currTB = m.during{tt,state}(m.XX,m.YY);
else
    currTB = m.duringRaw{state}(:,:,tt);
end
switch plotType
    case 'imagesc'
        imagesc(m.XX(:,1)',m.YY(1,:),currTB',ylims);
        set(gca,'YDir','normal');
    case 'surf'
        surf(m.XX,m.YY,currTB);view([45 45]);hold on
        XX = m.XX;
        YY = m.YY;
        points=[[min(XX(:)) min(XX(:)) max(XX(:)) max(XX(:))]' ...
            [0 max(YY(:)) max(YY(:)) 0]' (m.before{state}).*[1 1 1 1]'];
        h = fill3(points(:,1),points(:,2),points(:,3),'k');
        grid on
        alpha(h, 0.1)
        caxis([0 1]);
        zlim([0 1]);
end

colormap([0 0 0;jet]);
title([m.key{state} ' ' num2str(tGrid(tt)) ' s from FE'])
xlabel('dF (spks/s^2)');ylabel('F (spks/s)');

end

function [XAll,yAll,OLS,OLSdFOnly] = calcOLS(m,tt,state)
X = [ones(numel(m.XX),1) m.XX(:) m.YY(:)];
y = m.during{tt,state}(m.XX,m.YY);
y = y(:);
XAll = X;yAll = y;

OLS = inv(X'*X)*X'*y;
OLSdFOnly = inv(X(:,1:2)'*X(:,1:2))*X(:,1:2)'*y;
end

function [] = plotCB(m,ztick,ylims,state)

cb=colorbar('XTick', ztick);cb.Position = cb.Position + [0.6e-1, 0.1e-2, 0, 0];
cb.Ticks(cb.Ticks<ylims(1)) = [];
[newTicks,ndx] = sort([cb.Ticks, m.before{state}+eps]);
newLabels = [cb.TickLabels; 'baseline'];
cb.Ticks = newTicks;
cb.TickLabels = newLabels(ndx);

end








