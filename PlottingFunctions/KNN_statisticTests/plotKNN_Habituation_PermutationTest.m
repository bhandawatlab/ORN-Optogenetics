function [fNum,stateFigure,kinFigure] = plotKNN_Habituation_PermutationTest(f_orcoAll,state2Cons,tim2Cons,fNum)

k=ones(1,numel(state2Cons));
nGen_per_page = 3;
nPage_all_gen = ceil(numel(f_orcoAll)./nGen_per_page);
fNum_byState = fNum:nPage_all_gen:numel(state2Cons).*nPage_all_gen;
nGen = size(f_orcoAll,2);
for gen = 1:nGen
    m = f_orcoAll{gen}.model.params;
    tGrid = f_orcoAll{1,gen}.model.TurnBias.tt;
    XX = repmat(f_orcoAll{1,gen}.model.TurnBias.XX,1,1,size(m{1}.val,3));
    YY = repmat(f_orcoAll{1,gen}.model.TurnBias.YY,1,1,size(m{1}.val,3));
    
    for i = 1:numel(state2Cons)
        bootstrapped_h = m{i}.HabituationBootstrap.h;
        
        if mod(k(i),24) == 1
            figure(fNum_byState(i));set(gcf,'Position',[2 42 838 924])
            currSupT = [f_orcoAll{gen}.id];
            for g = 1:min(nGen-gen,2)
                currSupT = [currSupT '; ' f_orcoAll{gen+g}.id];
            end
            currSupT = {currSupT; 'Habituation Permutation Test 1=>,-1=< 95CI';' '};
            suptitle(currSupT)
            stateFigure{fNum_byState(i)} = m{i}.state.state;
            kinFigure{fNum_byState(i)} = m{i}.state.kin;
            
            k(i) = 1;fNum_byState(i) = fNum_byState(i)+1;
        end
        figure(fNum_byState(i)-1);
        
        for tt = 1:min(numel(tim2Cons),8)%24
            subplot(6,4,k(i));
            t = tim2Cons(tt)==tGrid;
            plotSubplot(XX, YY, bootstrapped_h, [-1.1 1], t);
            title([m{i}.state.state ' ' m{i}.state.kin ' ' num2str(tGrid(t)) ' s'])
            k(i)=k(i)+1;
        end
        
        if mod(k(i)-1,4)==0
            cb=colorbar('XTick', round([-1 0 1],2));
            cb.Position = cb.Position + [0.45e-1, 0.1e-2, 0, 0];
        end 
        if mod(k(1),4)~=1
            ylabel('')
        end
        k(i) = k(i)+(8-tt);
        
    end
    %k = k+(8-numel(state2Cons));
end

end

function [] = plotSubplot(XX, YY, bootstrapped_h, ylims, tt)

imagesc(XX(:,1)',YY(1,:),bootstrapped_h(:,:,tt)',[ylims(1) ylims(2)]);set(gca,'YDir','normal');
colormap([0 0 0;jet]);
%plotCB(m,ztick)
xlabel('dF (spks/s^2)');ylabel('F (spks/s)');
xlim([-150 150]);ylim([0 50]);
xticks([-150:75:150]);yticks([0:25:50]);

end