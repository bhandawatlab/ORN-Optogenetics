function [] = violinPlots(data,yl,plotdist,plotLegend)
% violinPlots  plots violin plots with box plots and data scatter
%
% Inputs:
%    data: vector of data points
%    yl: y-limit
%    plotdist: True/false of whether to plot the distribution
%    plotLegend: True/false of whether to plot the legend
%

if ~isempty(yl)
    if yl(1)>=0
        [f,~] = ksdensity(data,data,'Support','positive');
        [f_dens,x_dens] = ksdensity(data,linspace(yl(1),yl(2),200),'Support','positive');
    else
        [f,~] = ksdensity(data,data);
        [f_dens,x_dens] = ksdensity(data,linspace(yl(1),yl(2),200));
    end
else
    [f,~] = ksdensity(data,data);
    [f_dens,x_dens] = ksdensity(data);
end
m_f = max(f_dens);
f = f./m_f;
x = (rand(size(data))-0.5).*f;
scatter(x,data,'k','LineWidth',1);hold on;
plot([-0.5 0.5],mean(data)*[1 1],'-r')
plot([-0.5 0.5],median(data)*[1 1],'--r')
if plotdist
    if isempty(yl)
        yl = ylim;
    end
    
    plot([f_dens;-f_dens]./m_f./2,x_dens,'-b','linewidth',1)
    if plotLegend
        legend({'data','mean','median','Density'})
    end
else
    if plotLegend
        legend({'data','mean','median'})
    end
end
text(-0.5, yl(2).*0.9,['mean=' num2str(round(mean(data)))])
text(-0.5, yl(2).*0.8,['median=' num2str(round(median(data)))])
xlim([-0.5 0.5]);ylim(yl)
set(gca,'xtick',[])
hold off
end