function [fNum,stateFigure,kinFigure] = plotStatTestAdaptationKNN(f_orcoAll,state2Cons,tim2Cons,fNum)

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
        
        
        KNN_rawData = m{i}.KNN_rawData;
        nXY_grids = prod(size(XX,[1,2]));
        KNN_dataByTime = nan(nXY_grids,m{i}.KNN.K,numel(tGrid));
        for t_slice = 1:numel(tGrid)
            currNdx = (t_slice-1)*nXY_grids+1:t_slice*nXY_grids;
            KNN_dataByTime(:,:,t_slice)  = KNN_rawData(currNdx,:);
        end
        
        goodDat = sum(squeeze(sum(~isnan(KNN_dataByTime),[2]))>15,2);
        KS_test_h = nan(size(XX));
        KS_test_p = nan(size(XX));
        for t_slice = 2:numel(tGrid)
            KS_test_h_tmp = nan(size(KNN_dataByTime,1),1);
            KS_test_p_tmp = nan(size(KNN_dataByTime,1),1);
            goodDat = find(all(squeeze(sum(~isnan(KNN_dataByTime(:,:,[1 t_slice])),2))>15,2));
            for j = 1:numel(goodDat)
                tmp_dat1 = KNN_dataByTime(goodDat(j),:,1);
                tmp_dat1(isnan(tmp_dat1)) = [];
                tmp_dat2 = KNN_dataByTime(goodDat(j),:,t_slice);
                tmp_dat2(isnan(tmp_dat2)) = [];
                
                [h,p] = kstest2(tmp_dat1,tmp_dat2);
            
                KS_test_h_tmp(goodDat(j),1) = h;
                KS_test_p_tmp(goodDat(j),1) = p;
            end
            KS_test_h(:,:,t_slice) = reshape(KS_test_h_tmp,size(XX,[1,2]));
            KS_test_p(:,:,t_slice) = reshape(KS_test_p_tmp,size(XX,[1,2]));
        end
        
        log_p = log10(KS_test_p);
        log_p(log_p<-3) = -3;
        
        if mod(k(i),24) == 1
            figure(fNum_byState(i));set(gcf,'Position',[2 42 838 924])
            currSupT = [f_orcoAll{gen}.id];
            for g = 1:min(nGen-gen,2)
                currSupT = [currSupT '; ' f_orcoAll{gen+g}.id];
            end
            currSupT = [{currSupT; 'KS test from first entry (log10(prob))';' '}];
            suptitle(currSupT)
            stateFigure{fNum_byState(i)} = m{i}.state.state;
            kinFigure{fNum_byState(i)} = m{i}.state.kin;
            
            k(i) = 1;fNum_byState(i) = fNum_byState(i)+1;
        end
        figure(fNum_byState(i)-1);
        
        for tt = 1:min(numel(tim2Cons),8)%24
            subplot(6,4,k(i));
            t = tim2Cons(tt)==tGrid;
            plotSubplot(XX, YY, log_p, [-3.5 0], t);
            title([m{i}.state.state ' ' m{i}.state.kin ' ' num2str(tGrid(t)) ' s'])
            k(i)=k(i)+1;
        end
        
        if mod(k(i)-1,4)==0
            cb=colorbar('XTick', round(log10([0.001 0.01 0.05]),2));
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

function [] = plotSubplot(XX, YY, KS_test_p, ylims, tt)

imagesc(XX(:,1)',YY(1,:),KS_test_p(:,:,tt)',[ylims(1) ylims(2)]);set(gca,'YDir','normal');
colormap([0 0 0;jet]);
%plotCB(m,ztick)
xlabel('dF (spks/s^2)');ylabel('F (spks/s)');
xlim([-150 150]);ylim([0 50]);
xticks([-150:75:150]);yticks([0:25:50]);

end