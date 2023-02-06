function [fNum,stateFigure,kinFigure] = plotStatTestRetContKNN(f_orcoAll,state2Cons,tim2Cons,p_cutoff,fNum)

if isempty(p_cutoff)
    p_cutoff = nan;
end

k=1;
nGen_per_page = 3;
nGen = size(f_orcoAll,2);
for gen = 1:nGen
    m_1 = f_orcoAll{1,gen}.model.params;
    m_2 = f_orcoAll{2,gen}.model.params;
    tGrid = f_orcoAll{1,gen}.model.TurnBias.tt;
    XX = repmat(f_orcoAll{1,gen}.model.TurnBias.XX,1,1,size(m_1{1}.val,3));
    YY = repmat(f_orcoAll{1,gen}.model.TurnBias.YY,1,1,size(m_1{1}.val,3));
    
    for i = 1:numel(state2Cons)
        if mod(k,24) == 1
            figure(fNum);set(gcf,'Position',[2 42 838 924])
            currSupT = [f_orcoAll{1,gen}.id];
            currSupT_2 = [f_orcoAll{2,gen}.id];
            for g = 1:min(nGen-gen,2)
                currSupT = [currSupT '; ' f_orcoAll{1,gen+g}.id];
                currSupT_2 = [currSupT_2 '; ' f_orcoAll{2,gen+g}.id];
            end
            if isnan(p_cutoff)
                currSupT = {[currSupT ' vs '];[currSupT_2 ' KS test (log10(prob))'];' '};
            else
                currSupT = {[currSupT ' vs '];[currSupT_2 ' KS test p<' num2str(p_cutoff)];' '};
            end
            stateFigure{fNum} = 'all';
            kinFigure{fNum} = 'all';
            suptitle(currSupT)
            k(1) = 1;fNum = fNum+1;
        end
        
        KNN_rawData_1 = m_1{i}.KNN_rawData;
        KNN_rawData_2 = m_2{i}.KNN_rawData;
        
        try
            goodDat = find(sum(~isnan(KNN_rawData_1) & ~isnan(KNN_rawData_2),2)>15);
        catch
            [maxK,ndx] = max([size(KNN_rawData_1,2),size(KNN_rawData_2,2)]);
            if ndx == 2
                KNN_rawData_1 = [KNN_rawData_1 nan(size(KNN_rawData_1,1))];
                KNN_rawData_1 = KNN_rawData_1(:,1:maxK);
            else
                KNN_rawData_2 = [KNN_rawData_2 nan(size(KNN_rawData_2,1))];
                KNN_rawData_2 = KNN_rawData_2(:,1:maxK);
            end
            goodDat = find(sum(~isnan(KNN_rawData_1) & ~isnan(KNN_rawData_2),2)>15);
        end
        
        KS_test_h = nan(size(KNN_rawData_1,1),1);
        KS_test_p = nan(size(KNN_rawData_1,1),1);
        for j = 1:numel(goodDat)
            tmp_dat1 = KNN_rawData_1(goodDat(j),:);
            tmp_dat1(isnan(tmp_dat1)) = [];
            tmp_dat2 = KNN_rawData_2(goodDat(j),:);
            tmp_dat2(isnan(tmp_dat2)) = [];
            [h,p] = kstest2(tmp_dat1,tmp_dat2);
            
            KS_test_h(goodDat(j),1) = h;
            KS_test_p(goodDat(j),1) = p;
        end
        
        KS_test_h = reshape(KS_test_h,size(XX));
        KS_test_p = reshape(KS_test_p,size(XX));
        
        log_p = log10(KS_test_p);
        log_p(log_p<-3) = -3;
        
        tt = tim2Cons(1)==tGrid;
        subplot(6,4,k);
        if isnan(p_cutoff)
            plotSubplot(XX(:,:,1), YY(:,:,1), log_p, [-3.5 0], 273, tt);
        else
            KS_test_p_cutoff = double(KS_test_p<p_cutoff);
            KS_test_p_cutoff(isnan(KS_test_p)) = nan;
            plotSubplot(XX, YY, KS_test_p_cutoff, [-0.1 1], 150, tt);%205
        end
        title([m_1{i}.state.state ' ' m_1{i}.state.kin ' ' num2str(tGrid(tt)) ' s'])
        
        if mod(k(1),4)==0
            if isnan(p_cutoff)
                cb=colorbar('XTick', round(log10([0.001 0.01 0.05]),2));
            else
                cb=colorbar('XTick', [0 1]);
            end
            cb.Position = cb.Position + [0.45e-1, 0.1e-2, 0, 0];
        end 
        if mod(k(1),4)~=1
            ylabel('')
        end
        k(1) = k(1)+1;
    end
    k = k+(8-numel(state2Cons));
end

end

function [] = plotSubplot(XX, YY, KS_test_p, ylims, n, tt)
cc = jet(272);

imagesc(XX(:,1)',YY(1,:),KS_test_p(:,:,tt)',[ylims(1) ylims(2)]);set(gca,'YDir','normal');
colormap([0 0 0; cc(end-n:end-10,:)]);
%plotCB(m,ztick)
xlabel('dF (spks/s^2)');ylabel('F (spks/s)');
xlim([-150 150]);ylim([0 50]);
xticks([-150:75:150]);yticks([0:25:50]);

end