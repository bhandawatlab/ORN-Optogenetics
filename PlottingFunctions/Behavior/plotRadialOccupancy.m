function [fNum] = plotRadialOccupancy(f_orcoAll,lab,border,stopSpd,dur,fNum)

yl = 0.4;
k = 1;
for gen = 1:numel(f_orcoAll)
    f_orco = f_orcoAll{gen};
    
    if ~isempty(f_orco)
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
        plot(radialProb.before.x.*f_orco.rBound,radialProb.before.y,'-g','Linewidth',1);hold on;
        plot(radialProb.during.x.*f_orco.rBound,radialProb.during.y,'-r','Linewidth',1);
        plot([border,border],[0 yl],'--k');hold off;
        text(0.05*f_orco.rBound,0.3,['n=' num2str(f_orco.nFly)])
        ylim([0 yl]);xlim([0 f_orco.rBound]);
        xlabel('radial pos (cm)');ylabel('Probability');
        xticks([0:0.25:1].*f_orco.rBound);yticks([0:0.1:yl]);
        
        if k==1
            legend({'Before','During'})
        end
        
        title(f_orco.id);
    end
    k = k+1;
end

end