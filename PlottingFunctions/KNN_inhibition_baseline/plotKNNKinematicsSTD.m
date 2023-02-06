function [fNum,stateFigure,kinFigure] = plotKNNKinematicsSTD(f_orcoAll,state2Cons,tim2Cons,plotType,colortype,fNum)

% loop through each genotype
k=ones(1,numel(state2Cons));
nGen_per_page = 3;
nPage_all_gen = ceil(numel(f_orcoAll)./nGen_per_page);
fNum_byState = fNum:nPage_all_gen:numel(state2Cons).*nPage_all_gen;
nGen = numel(f_orcoAll);
for gen = 1:nGen
    f_orco = f_orcoAll{gen};
    m = f_orco.model.TurnBias;
    tGrid = m.tt;
    XX = f_orco.model.TurnBias.XX;
    YY = f_orco.model.TurnBias.YY;
    
    for i = 1:numel(state2Cons)
        state = state2Cons(i);
        m = f_orco.model.params{state};
        
        for slice = 1:numel(m.g2_std)
            m.val(:,:,slice) = m.g2_std{slice}(XX,YY);
        end
        
        if strcmpi(m.state.kin,'spd')
            ylab = 'mm/s';ylims = [0.5 1.5];ztick = [0.5 1 1.5];
        elseif strcmpi(m.state.kin,'totCurv')
            ylab = 'deg';ylims = [0 2];ztick = [0 1 2];
        elseif strcmpi(m.state.kin,'avgCurv')
            ylab = 'deg/s';ylims = [0 3];ztick = [0 1.25 2.5];
        elseif strcmpi(m.state.kin,'dur') && strcmpi(m.state.state,'stops')
            ylab = 's';ylims = [0 5];ztick = [0 2.5 5];
        elseif strcmpi(m.state.kin,'dur') && ~strcmpi(m.state.state,'stops')
            ylab = 's';ylims = [0.5 1.5];ztick = [0.5 1 1.5];
        end
        
        % average across time
        if numel(tim2Cons)==1
            if mod(k(1),24) == 1
                figure(fNum);set(gcf,'Position',[2 42 838 924])
                currSupT = [f_orcoAll{gen}.id];
                for g = 1:min(nGen-gen,2)
                    currSupT = [currSupT '; ' f_orcoAll{gen+g}.id];
                end
                currSupT = [currSupT '; time from first entry'];
                suptitle(currSupT)
                stateFigure{fNum} = 'all';
                kinFigure{fNum} = 'all';
                k(1) = 1;fNum = fNum+1;
            end
            
            subplot(6,4,k(1));
            plotSubplot(XX, YY, m, ylims, tim2Cons(1)==tGrid, tGrid, ...
                ztick, colortype, plotType);
            plotCB(m,ztick)
            if mod(k(1),4)~=1
                ylabel('')
            end
            k(1) = k(1)+1;
        else% across time
            if mod(k(state),24) == 1
                figure(fNum_byState(state));set(gcf,'Position',[2 42 838 924])
                currSupT = [f_orcoAll{gen}.id];
                for g = 1:min(nGen-gen,2)
                    currSupT = [currSupT '; ' f_orcoAll{gen+g}.id];
                end
                currSupT = [currSupT '; time from first entry'];
                suptitle(currSupT)
                stateFigure{fNum_byState(state)} = m.state.state;
                kinFigure{fNum_byState(state)} = m.state.kin;
                
                k(state) = 1;fNum_byState(state) = fNum_byState(state)+1;
            end
            figure(fNum_byState(state)-1);
            
            for tt = 1:min(numel(tim2Cons),8)%24
                subplot(6,4,k(state));
                plotSubplot(XX, YY, m, ylims, tim2Cons(tt)==tGrid, tGrid, ...
                    ztick, colortype, plotType);
                k(state)=k(state)+1;
            end
            k(state) = k(state)+(8-tt);
            
            if mod(k(state)-1,8)==0
                plotCB(m,ztick)
            end
        end
    end
    
    if numel(tim2Cons)==1
        k(1) = k(1)+(8-numel(state2Cons));
    end
end

fNum = fNum_byState(end);

end

function [] = plotSubplot(XX, YY, m, ylims, tt, tGrid, ztick, colortype, plotType)
switch colortype
    case 'absolute'
        ylims(1) = ylims(1);
    case 'relative'
        baseline = (m.baseline(2));
        if baseline<=mean(ylims)
            ylims(1) = baseline-(ylims(2)-baseline);
        else
            ylims(2) = baseline.*2;
        end
end

switch plotType
    case 'imagesc'
        imagesc(XX(:,1)',YY(1,:),(m.val(:,:,tt))',[ylims(1) ylims(2)]);set(gca,'YDir','normal');
        colormap([0 0 0;jet]);
        %plotCB(m,ztick)
    case 'surf'
        surf(XX,YY,(m.val(:,:,tt)));view([45 45]);hold on;
        points=[[min(XX(:)) min(XX(:)) max(XX(:)) max(XX(:))]' ...
            [0 max(YY(:)) max(YY(:)) 0]' (m.baseline(2)).*[1 1 1 1]'];
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

title([m.state.state ' ' m.state.kin ' ' num2str(tGrid(tt)) ' s'])
xlabel('dF (spks/s^2)');ylabel('F (spks/s)');
xlim([-150 150]);ylim([0 50]);
xticks([-150:75:150]);yticks([0:25:50]);

end

function [] = plotCB(m,ztick)

cb=colorbar('XTick', ztick);cb.Position = cb.Position + [0.45e-1, 0.1e-2, 0, 0];
[newTicks,ndx] = sort([cb.Ticks, (m.baseline(2))]);
cb.TickLabels = cellstr(num2str(cb.Ticks'));
newLabels = [cb.TickLabels; 'baseline'];
cb.Ticks = newTicks;
cb.TickLabels = newLabels(ndx);

end