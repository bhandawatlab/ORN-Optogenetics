function [fNum] = plotDistributionComparisons(f_orcoAll,state2Cons,fNum)

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
    
    for state = 1:numel(state2Cons)
        for i = 1:numel(distTypes)
            for j = 1:numel(distTypes)
                % note that likelihood ratio test is not valid because
                % these distributions are not nested
                %L_ratio_test{state,gen}(i,j) = -2.*(m{state}.KNN.nLL(i)-m{state}.KNN.nLL(j));
                
                % -log(L(null)/L(test)) = nLL(null)-nLL(test)
                % positive means test works better
                L_ratio{state,gen}(i,j) = (m{state}.KNN.nLL(i)-m{state}.KNN.nLL(j));
            end
        end
    end
end

for i = 1:numel(distTypes)
    for j = (i+1):numel(distTypes)
        figure(fNum);set(gcf,'Position',[2 42 838 924])
        k = 1;
        for state = 1:numel(state2Cons)
            subplot(3,3,k);
            tmp = cellfun(@(x) x(i,j),L_ratio(state,:));
            %histogram(tmp,3);
            if k~= numel(state2Cons)
                violinPlotsStats(tmp,[]);
                %violinPlots(tmp,[],true,false);
            else
                violinPlots(tmp,[],true,true);
            end
            hold on;plot([-0.5 0.5],[0 0],'k');hold off
            title([m{state}.state.state ' ' m{state}.state.kin])
            k = k+1;
        end
        sgtitle(['nLL(' distTypes{i} ') - nLL(' distTypes{j} ')'])
        fNum = fNum+1;
    end
end

% 
% for state = 1:numel(state2Cons)
%     figure(fNum);set(gcf,'Position',[2 42 838 924])
%     k = 1;
%     for i = 1:numel(distTypes)
%         for j = (i+1):numel(distTypes)
%             subplot(4,3,k);
%             tmp = cellfun(@(x) x(i,j),L_ratio(state,:));
%             %histogram(tmp,3);
%             if k~= numel(distTypes)*(numel(distTypes)-1)/2
%                 violinPlots(tmp,[],true,false);
%             else
%                 violinPlots(tmp,[],true,true);
%             end
%             hold on;plot([-0.5 0.5],[0 0],'k');hold off
%             title({['null=' distTypes{i}],['test=' distTypes{j}]})
%             k = k+1;
%         end
%     end
%     sgtitle([m{state}.state.state ' ' m{state}.state.kin])
%     fNum = fNum+1;
% end

end