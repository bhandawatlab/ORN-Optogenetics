function [] = plotSignificanceChangeByROI(allMat,type,thresh,plotSyn)

SynergyIdx = 11;
MuIdx = [4:5];

%plot
nGen2Cons = 3;

pBaseline = zeros(nGen2Cons);
pGen(1,:)  = [1,0];
pGen(2,:)  = [0,1];
pGen(3,:)  = [1,1];
pGen = [pBaseline,pGen];
pCovGen(:,1) = [0;0;1];

x = [2,8,4,6,5];

dFType = {'-df','0df','+df'};
fType = {'highF','lowF','Inh'};

states2Cons = [1,2,3,4];
for i = 1:numel(allMat)
    figure;set(gcf,'Position',[2 42 838 924])
    load(allMat{i},'allMedBestFit','f_orco','RegionLabel','labCov')
    if strcmpi(type,'absolute')
        colLab = {['mu ' labCov{1}],['mu ' labCov{2}],'muAB','absolute interactions'};
    elseif strcmpi(type,'relative')
        colLab = {['mu ' labCov{1}],['mu ' labCov{2}],'muAB','relative interactions'};
    end
    covar = zeros(numel(states2Cons),5);
    
    varAB = cell(numel(states2Cons),4);
    synergy = cell(numel(states2Cons),4);
    k = 1;
    for state = states2Cons
        varAB(state,:) = {zeros(3,3)};
        synergy(state,:) = {zeros(3,3)};
        m = f_orco.model.params{state};
        stateType{1,state} = m.state.state;
        stateType{2,state} = m.state.kin;
        for ROI = 1:5
            tmp = exp(allMedBestFit{state,ROI}(1:5)*pGen'+allMedBestFit{state,ROI}(11)*pCovGen');
            for j = 1:3
                varAB{state,j}(x(ROI)) = tmp(j);
            end
            if strcmpi(type,'absolute')
                synergy{state,4}(x(ROI)) = tmp(3)*(1-1./exp(allMedBestFit{state,ROI}(SynergyIdx)));
            elseif strcmpi(type,'relative')
                synergy{state,4}(x(ROI)) = exp(allMedBestFit{state,ROI}(SynergyIdx));
            end
            if plotSyn
                prodMat = sign(prod(allMedBestFit{state,ROI}(MuIdx)));
                synMat = sign((allMedBestFit{state,ROI}(SynergyIdx)));
                singleMat = sign(allMedBestFit{state,ROI}(MuIdx(1)));
                
                if strcmpi(type,'absolute')
                    sigChange = abs(synergy{state,4}(x(ROI)))>thresh(state);
                elseif strcmpi(type,'relative')
                    sigChange = abs(synergy{state,4}(x(ROI))-1)>thresh(state);
                end
                
                if prodMat==1 && synMat==singleMat && sigChange
                    varAB{state,4}(x(ROI)) = 1;
                elseif prodMat==1 && synMat~=singleMat && sigChange
                    varAB{state,4}(x(ROI)) = -1;
                elseif ~sigChange
                    varAB{state,4}(x(ROI)) = -1.5;
                else
                    varAB{state,4}(x(ROI)) = 0.01;
                end
            end
        end
        
        for j = 1:4
            subplot(numel(states2Cons),4,k)
            tmp =varAB{state,j};
            tmp(setdiff(1:9,x)) = -1.5;
            imagesc(tmp,[-1.5 1])
            
            set(gca,'ytick',[1:1:8],'yticklabel',fType);
            set(gca,'xtick',[1:1:5],'xticklabel',dFType);
            cMap = jet(221);
            %cMap(11,:)=1;%101
            cMap=[ones(50,3);cMap];%101
            colormap(cMap);
            
            if j <4
                t2 = compose('%g',round(varAB{state,j},2));
            else
                t2 = compose('%g',round(synergy{state,j},2));
            end
            
            [xx,yy] = meshgrid([1:3],[1:3]);
            text(xx(:), yy(:), t2(:), 'HorizontalAlignment', 'Center','Color',0.5.*[1 1 1])
            title({[stateType{1,state} ' ' stateType{2,state}];colLab{j}})
            k = k+1;
        end
        cb=colorbar;cb.Position = cb.Position + [0.75e-1, 0, 0, 0];
        cb.YTick = [-1.5 -1 0 1];
        cb.YTickLabel = {'No Interactions' ,'Ant', 'N/A', 'Syn'};
    end
    sgtitle([labCov{1,1} ' + ' labCov{2,1} ' in real space (mm/s,deg,deg/s)'])
end

end
