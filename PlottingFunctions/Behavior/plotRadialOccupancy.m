function [fNum] = plotRadialOccupancy(f_orcoAll,lab,border,stopSpd,dur,fNum)

yl = 0.4;
k = 1;
for gen = 1:numel(f_orcoAll)
    f_orco = f_orcoAll{gen};
    
    if ~isempty(f_orco)
        if gen>1
            radialProbPrevious = radialProb;
        end
        radialProb = f_orco.getRadialProb(lab,border,stopSpd,dur);
        
        if k ==1 || k > 24
            figure(fNum);set(gcf,'Position',[2 42 838 924])
            if strcmpi(lab,'H')
                suptitle('Radial Occupancy: Head')
            else
                suptitle('Radial Occupancy: Head')
            end
            k = 1;fNum = fNum+1;
        end
        
        subplot(6,4,k); 
        p_b = shadedErrorBar(radialProb.before.x.*f_orco.rBound,radialProb.before.weightedMu,...
            radialProb.before.sem,'lineprops',{'-g','Linewidth',1});hold on;
        p_d = shadedErrorBar(radialProb.during.x.*f_orco.rBound,radialProb.during.weightedMu,...
            radialProb.during.sem,'lineprops',{'-r','Linewidth',1});
        [~,p] = kstest2(cell2mat(radialProb.before.yRaw),cell2mat(radialProb.during.yRaw));
        text(0.05,0.35,['ks bef dur p=' num2str(p)], 'Interpreter', 'none')
        if mod(gen,2)==0
            [~,p] = kstest2(cell2mat(radialProb.during.yRaw),cell2mat(radialProbPrevious.during.yRaw));
            text(0.05,0.3,['ks cont ret p=' num2str(p)], 'Interpreter', 'none')
        end
%         plot(radialProb.before.x.*f_orco.rBound,radialProb.before.y,'-g','Linewidth',1);hold on;
%         plot(radialProb.during.x.*f_orco.rBound,radialProb.during.y,'-r','Linewidth',1);
        p_border = plot([border,border],[0 yl],'--k');hold off;
        text(0.05*f_orco.rBound,0.25,['n=' num2str(f_orco.nFly)])
        ylim([0 yl]);xlim([0 f_orco.rBound]);
        xlabel('radial pos (cm)');ylabel('Probability');
        xticks([0:0.25:1].*f_orco.rBound);yticks([0:0.1:yl]);
        
        if k==1
            legend([p_b.mainLine p_d.mainLine p_border],{'Before','During','border'})
        end
        
        title(f_orco.id);
    end
    k = k+1;
end

end