function [fNum,stateFigure,kinFigure] = plotDistributionComparisons2(f_orcoAll,state2Cons,fNum)

k=1;
nGen_per_page = 3;
nPage_all_gen = ceil(numel(f_orcoAll)./nGen_per_page);
nGen = numel(f_orcoAll);
for gen = 1:nGen
    f_orco = f_orcoAll{gen};
    m = f_orco.model.params;
    
    distTypes = m{1}.KNN.DistTypes;
    distTypes = strrep(distTypes,'Normal', 'Norm');
    distTypes = strrep(distTypes,'Lognormal', 'logN');
    distTypes = strrep(distTypes,'Gama', 'gam');
    distTypes = strrep(distTypes,'Weibull', 'wb');
    distTypes = strrep(distTypes,'Exponential', 'exp');
    
    
    for i = 1:numel(state2Cons)
        if mod(k,24) == 1
            figure(fNum);set(gcf,'Position',[2 42 838 924])
            currSupT = [f_orcoAll{gen}.id];
            for g = 1:min(nGen-gen,2)
                currSupT = [currSupT '; ' f_orcoAll{gen+g}.id];
            end
            currSupT = [currSupT '; Distribution Fits'];
            stateFigure{fNum} = 'all';
            kinFigure{fNum} = 'all';
            suptitle(currSupT)
            k(1) = 1;fNum = fNum+1;
        end
        
        subplot(6,4,k);
        yyaxis left; plot(m{i}.KNN.nLL,'-k','LineWidth',1);hold on;
        yyaxis right; plot(m{i}.KNN.AIC,'-r','LineWidth',1);
        plot(m{i}.KNN.BIC,'-b','LineWidth',1);
        
        xticks(1:length(m{i}.KNN.DistTypes))
        xticklabels(distTypes)
        xtickangle(45)
        title([m{i}.state.state ' ' m{i}.state.kin])
        
        if mod(k,24) == 1
            legend({'-LL','AIC','BIC'})
        end
        if mod(k,4)==1
            yyaxis left; ylabel('neg LL')
        elseif mod(k,4)==0
            yyaxis right; ylabel('AIC or BIC')
        end
        k = k+1;
    end
    k = k+(8-numel(state2Cons));
end



end