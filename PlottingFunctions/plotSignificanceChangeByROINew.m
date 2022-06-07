function [allData] = plotSignificanceChangeByROINew(allMat,t1,t2,t3,baseline,type,thresh)

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


ROIPos{1} = [-200,0,t1+200,50];%negative df
ROIPos{2} = [t2,0,200-t2,50];%positive df
ROIPos{3} = [t1,t3,t2-t1,50-t3];%high F
ROIPos{4} = [t1,0,t2-t1,baseline];%inhibition
ROIPos{5} = [t1,baseline,t2-t1,t3];%low F
ROILabel = {'-df','+df','high F', 'inh','low F'};

%cc = [0.8500 0.3250 0.0980;0.3010 0.7450 0.9330;1 1 1;[213 196 161]./255];

cc = [0.8500 0.3250 0.0980;0.3010 0.7450 0.9330;1 1 1;[213 196 161]./255;...
    0.4660 0.6740 0.1880;0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];

dFType = {'-df','0df','+df'};
fType = {'highF','lowF','Inh'};

allData = cell(1,numel(allMat));

states2Cons = [1,2,3,4];
for i = 1:numel(allMat)
    figure;set(gcf,'Position',[2 42 838 924])
    load(allMat{i},'allMedBestFit','f_orco','RegionLabel','labCov')
    if strcmpi(type,'absolute')
        colLab = {['mu ' labCov{1}],['mu ' labCov{2}],'muAB','absolute interactions','dominate'};
    elseif strcmpi(type,'relative')
        colLab = {['mu ' labCov{1}],['mu ' labCov{2}],'muAB','relative interactions','dominate'};
    end    
    values = cell(numel(states2Cons),4);
    contributions = cell(numel(states2Cons),4);
    k = 1;
    for state = states2Cons
        values(state,:) = {zeros(size(RegionLabel))};
        contributions(state,:) = {zeros(size(RegionLabel))};
        m = f_orco.model.params{state};
        stateType{1,state} = m.state.state;
        stateType{2,state} = m.state.kin;
        for ROI = 1:5
            tmp = exp(allMedBestFit{state,ROI}(1:5)*pGen'+allMedBestFit{state,ROI}(11)*pCovGen');
            % genotype term
            for j = 1:3
                values{state,j}(ROI) = tmp(j);
            end
            % interaction term
            if strcmpi(type,'absolute')
                values{state,4}(ROI) = tmp(3)*(1-1./exp(allMedBestFit{state,ROI}(SynergyIdx)));
            elseif strcmpi(type,'relative')
                values{state,4}(ROI) = exp(allMedBestFit{state,ROI}(SynergyIdx));
            end
            
            % type of contributions of each genotype
            for j = 1:4
                % check if the change is significant
                if strcmpi(type,'absolute')
                    sigChange = abs(values{state,j}(ROI))>thresh(state);
                    posChange = values{state,j}(ROI)>thresh(state);
                elseif strcmpi(type,'relative')
                    sigChange = abs(values{state,j}(ROI)-1)>thresh(state);
                    posChange = (values{state,j}(ROI)-1)>thresh(state);
                end
                if sigChange && posChange
                    % increase
                    contributions{state,j}(ROI) = 1;
                elseif sigChange && ~posChange
                    % decrease
                    contributions{state,j}(ROI) = 2;
                else
                    % no change
                    contributions{state,j}(ROI) = 3;
                end
            end
            
            
            % type of contributions of interactions
            % -------------------------------------------------------------
            % these values help to determine is something is synergistic or
            % antagonistic.
            prodMat = sign(prod(allMedBestFit{state,ROI}(MuIdx)));
            synMat = sign((allMedBestFit{state,ROI}(SynergyIdx)));
            singleMat = sign(allMedBestFit{state,ROI}(MuIdx(1)));
            
            if prodMat==1 && synMat==singleMat && sigChange
                % synergistic
                contributions{state,4}(ROI) = 1;
            elseif prodMat==1 && synMat~=singleMat && sigChange
                % antagonistic
                contributions{state,4}(ROI) = 2;
            elseif ~sigChange
                % no change
                contributions{state,4}(ROI) = 3;
            else
                % interaction, but not synergistic or antagonistic
                contributions{state,4}(ROI) = 4;
            end
            % -------------------------------------------------------------
            
%             % type of domination
%             % -------------------------------------------------------------
%             if contributions{state,4}(ROI)>2
%                 diff_AtoAB = abs(values{state,3}(ROI)-values{state,1}(ROI));
%                 diff_BtoAB = abs(values{state,3}(ROI)-values{state,2}(ROI));
%                 diff_AtBtoAB = abs(values{state,3}(ROI)-values{state,1}(ROI).*values{state,2}(ROI));
%                 diff_InttoAB = abs(values{state,3}(ROI)-values{state,4}(ROI));
%                 A_sig = contributions{state,1}(ROI)<3;
%                 B_sig = contributions{state,2}(ROI)<3;
%                 AB_sig = contributions{state,2}(ROI)<3;
%                 
%                 contributions{state,j}
%                 
%             end
%             % -------------------------------------------------------------
            
            % type of domination
            % -------------------------------------------------------------
            if contributions{state,4}(ROI)>=3
                A_B_opp = contributions{state,1}(ROI)==1 && contributions{state,2}(ROI)==2;
                A_B_opp_2 = contributions{state,1}(ROI)==2 && contributions{state,2}(ROI)==1;
                A_B_opp = A_B_opp | A_B_opp_2;
                if contributions{state,4}(ROI)==3 && A_B_opp
                    % if no interaction, and effect of A and B are
                    % different, then it's a linear mixture
                    contributions{state,4}(ROI) = 3;
                else
                    if (contributions{state,1}(ROI)~=contributions{state,2}(ROI)) %&& contributions{state,3}(ROI)~=3
                        if contributions{state,1}(ROI) == contributions{state,3}(ROI)
                            % dominated by A
                            contributions{state,4}(ROI) = 5;
                        elseif contributions{state,2}(ROI) == contributions{state,3}(ROI)
                            % dominated by B
                            contributions{state,4}(ROI) = 6;
                        end
    %                 else
    %                     % no domination
    %                     contributions{state,4}(ROI) = 7;
                    end
                end
                % -------------------------------------------------------------
            end
        end
        
        for j = 1:4
            subplot(numel(states2Cons),4,k)
            for ROI = 1:5
                rectangle('Position',ROIPos{ROI},'FaceColor',...
                    cc(contributions{state,j}(ROI),:),'EdgeColor','k',...
                    'LineWidth',1)
                txtLoc = ROIPos{ROI}(1:2)+ROIPos{ROI}(3:4)./2;
                if j<=4
                    t2 = num2str(round(values{state,j}(ROI),2));
                    text(txtLoc(1),txtLoc(2),t2, 'HorizontalAlignment', 'Center','Color',[0 0 0])
                end
            end
            xlim([-200 200]);ylim([0 50])
            set(gca,'ytick',[0, 10, 35],'yticklabel',fType);
            set(gca,'xtick',[-100, 0, 100],'xticklabel',dFType);
            
            title({[stateType{1,state} ' ' stateType{2,state}];colLab{j}})
            k = k+1;
        end
        colormap(cc)
        cb=colorbar;cb.Position = cb.Position + [0.5e-1, 0, 0, 0];
        cb.YTick = [0:(1./7):1]+(1./14);
        cb.YTickLabel = {'+/Syn' ,'-/Ant', 'no inter./lin. integate', 'other int',...
            [labCov{1,1} ' dom'],[labCov{2,1} ' dom']};
    end
    sgtitle([labCov{1,1} ' + ' labCov{2,1} ' in real space (mm/s,deg,deg/s)'])
    
    allData{i}.values = values;
    allData{i}.contributions = contributions;
    allData{i}.labCov = labCov;
    allData{i}.ROI = ROILabel;
    allData{i}.stateType = stateType;
end

end
