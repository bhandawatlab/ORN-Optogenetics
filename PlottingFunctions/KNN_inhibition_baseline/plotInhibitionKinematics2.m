function [fNum,stateFigure,kinFigure] = plotInhibitionKinematics2(f_orcoAll,state2Cons,tim2Cons,fNum)

t = 0:0.1:2;
k=ones(1,numel(state2Cons));
nGen_per_page = 3;
nPage_all_gen = ceil(numel(f_orcoAll)./nGen_per_page);
fNum_byState = fNum:nPage_all_gen:numel(state2Cons).*nPage_all_gen;
nGen = numel(f_orcoAll);

for gen = 1:nGen
    f_orco = f_orcoAll{gen};
    % plotting kinematics during inhibition
    
    m = f_orco.model.TurnBias;
    tGrid = m.tt;
    for i = 1:numel(state2Cons)
        state = state2Cons(i);
        m = f_orco.model.params{state};
        %figure(fNum);set(gcf,'Position',[2 42 838 924])
        %suptitle([f_orco.id ' ' m.state.state ' ' m.state.kin ' time from first entry'])
        if strcmpi(m.state.kin,'spd')
            ylab = 'mm/s';ylims = [0 20];ytick = [0 10 20];
        elseif strcmpi(m.state.kin,'totCurv')
            ylab = 'deg';ylims = [0 180];ytick = [0 90 180];
        elseif strcmpi(m.state.kin,'avgCurv')
            ylab = 'deg/s';ylims = [0 90];ytick = [0:15:45];
        end
        
        
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
            ttt = tim2Cons(1)==tGrid;
            
            subplot(6,4,k(1));
            plotSubplot(m,ttt,t,tGrid,ylims,ylab)
            k(1) = k(1)+1;
            
        else
            if mod(k(i),24) == 1
                figure(fNum_byState(i));set(gcf,'Position',[2 42 838 924])
                currSupT = [f_orcoAll{gen}.id];
                for g = 1:min(nGen-gen,2)
                    currSupT = [currSupT '; ' f_orcoAll{gen+g}.id];
                end
                currSupT = [currSupT '; time from first entry'];
                suptitle(currSupT)
                
                stateFigure{fNum_byState(i)} = m.state.state;
                kinFigure{fNum_byState(i)} = m.state.kin;

                k(i) = 1;fNum_byState(i) = fNum_byState(i)+1;
            end
            figure(fNum_byState(i)-1);

            for tt = 1:min(numel(tim2Cons),8)%24
                subplot(6,4,k(i));
                ttt = tim2Cons(tt)==tGrid;
                plotSubplot(m,ttt,t,tGrid,ylims,ylab)
                k(i)=k(i)+1;
            end
            k(state) = k(state)+(8-tt);
        end
    end
    if numel(tim2Cons)==1
        k(1) = k(1)+(8-numel(state2Cons));
    end
end

fNum = fNum_byState(end);
end

function [] = plotSubplot(m,tt,t,tGrid,ylims,ylab)

if ~isempty(m.inhibitionKin.a{tt})
    a = m.inhibitionKin.a{tt}(t);
    b = m.inhibitionKin.b{tt}(t);
    sampInh = [];
    if strcmpi(m.inhibitionKin.fitFun,'Lognormal') || strcmpi(m.inhibitionKin.fitFun,'logN')
        %--------------------log normal
        meanInh = exp(a);
        stdInh = exp([(a-b);(a+b)]);
    elseif strcmpi(m.inhibitionKin.fitFun,'beta')
        %--------------------beta
        meanInh = a./(a+b).*m.inhibitionKin.maxVal;
        stdInh = sqrt(a.*b./(((a+b).^2).*(a+b+1))).*m.inhibitionKin.maxVal;
        stdInh = meanInh+[-stdInh;stdInh];
    end
else
    % use KNN space, if not enough data for inhibition
    meanInh = exp(m.g2{tt}(0,0));
    stdInh = exp(m.g2{tt}(0,0)+[-(m.g2_std{tt}(0,0));(m.g2_std{tt}(0,0))]);
    meanInh = repmat(meanInh,size(t));
    stdInh = repmat(stdInh,size(t));
    text(0.1,ylims(2).*0.8,{'not enough data', 'using KNN'});hold on;
end
plot(t,meanInh,'k','Linewidth',2);hold on;
plot(t,stdInh,'--k','Linewidth',1);hold off
title([m.state.state ' ' m.state.kin ' ' num2str(tGrid(tt)) ' s'])
xlabel('time since inh (s)');ylabel(ylab)
xlim([0 max(t)]);ylim(ylims)
if tt == 1
    legend({'mu','mu-sig','mu+sig'})
end
end