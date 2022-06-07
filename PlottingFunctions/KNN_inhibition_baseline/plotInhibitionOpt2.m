function [fNum] = plotInhibitionOpt2(f_orcoAll,state,tim2Cons,fNum)

k = 1;
nGen = numel(f_orcoAll);
t0 = zeros(nGen,min(numel(tim2Cons),24));
for gen = 1:nGen
    f_orco = f_orcoAll{gen};
    % for border choice
    t = 0:0.1:2;
    m = f_orco.model.TurnBias;
    tGrid = m.tt;
    
    %figure(fNum);set(gcf,'Position',[2 42 838 924])
    for tt = 1:min(numel(tim2Cons),8)
        if mod(k,24) == 1
            figure(fNum);set(gcf,'Position',[2 42 838 924])
            currSupT = [f_orcoAll{gen}.id];
            for g = 1:min(nGen-gen,2)
                currSupT = [currSupT '; ' f_orcoAll{gen+g}.id];
            end
            currSupT = [currSupT ' Turn Optimality time from first entry'];
            suptitle(currSupT)
            k = 1;fNum = fNum+1;
        end
        tSlice = tGrid==tim2Cons(tt);
        
        subplot(6,4,k);
        if ~isempty(m.inhibitionKin{state}{tSlice})
            plot(t,m.inhibitionKin{state}{tSlice}(t),'k','Linewidth',1);
        else
            % use KNN space, if not enough data for inhibition
            KNN_inh = m.during{tSlice,state}(0,0);
            plot(t,repmat(KNN_inh,1,numel(t)),'k','Linewidth',1);
            text(0.1,0.2,{'not enough data', 'using KNN'})
        end
        title([num2str(tGrid(tSlice)) ' s'])
        xlabel('time since inh (s)');ylabel('Turn Optimality')
        xlim([0 max(t)]);ylim([0 1])
        
        k = k+1;
    end
    if numel(tim2Cons)<8
        k = k+(8-numel(tim2Cons));
    end
end

end