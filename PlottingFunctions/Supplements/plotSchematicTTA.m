function [] = plotSchematicTTA(f_orco)
%figureFile = 'TurnTriggeredAverageExampleTracks';

% get the sharp turn and stop indexes
stNdx = find(strcmpi(f_orco.states.key,'sharp turns'));
stopNdx = find(strcmpi(f_orco.states.key,'stops'));
% sample timeframe
tt = 5550:6435;

% example flies (note that fly 1 was used in Figure 2S4)
for fly = 1:8
    % find the start and endx frame for each sharp turn trajectory
    [startNdx_ST,endNdx_ST,type] = startEndSeq(f_orco.states.ndx(fly,tt)==stNdx);
    startNdx_ST(type==0) = [];
    endNdx_ST(type==0) = [];
    % find the start and endx frame for each stop trajectory
    [startNdx_stop,endNdx_stop,type] = startEndSeq(f_orco.states.ndx(fly,tt)==stopNdx);
    startNdx_stop(type==0) = [];
    endNdx_stop(type==0) = [];
    
    % index the x,y body and head position, firing rate, and curvature for
    % the period of time defined by tt
    currX = f_orco.x(fly,tt);
    currY = f_orco.y(fly,tt);
    currXH = f_orco.xH(fly,tt);
    currYH = f_orco.yH(fly,tt);
    currF = f_orco.spk(fly,tt);
    curv = f_orco.curv(fly,tt).*180./pi.*f_orco.fs;
    curv(f_orco.states.ndx(fly,tt)==stopNdx) = 0;%set the curvature at stops to be 0
    
    figure;set(gcf,'Position',[2 42 838 924]);
    subplot(4,1,[1 2]);hold on;
    % plot the body vectors
    for i = 1:numel(currY)
        plot([currX(i) currXH(i)],[currY(i) currYH(i)],'color',[0.5 0.5 0.5]);
    end
    % plot the x,y path trajectory
    %plot3(currX,currY,-tt./f_orco.fs,'k','LineWidth',1);view(2);hold on;
    plot(currX,currY,'k','LineWidth',1);
    % plot the sharp turn instances
    for i = 1:numel(startNdx_ST)
        plot(currX(startNdx_ST(i):endNdx_ST(i)),currY(startNdx_ST(i):endNdx_ST(i)),'r','LineWidth',2);
    end
    xlabel('x pos (cm)');ylabel('y pos (cm)')
    plotCircle([0 0],1.25,100,'c');
    plotCircle([0 0],4,100,'k');
    axis equal;
    
    % plot the stop instances
    for i = 1:numel(startNdx_stop)
        plot(currX(startNdx_stop(i):endNdx_stop(i)),currY(startNdx_stop(i):endNdx_stop(i)),'g','LineWidth',2);
    end
    
    % plot the curvature of the trajectory
    subplot(4,1,3);
    plot(tt./f_orco.fs,curv,'k','LineWidth',1);hold on;
    for i = 1:numel(startNdx_ST)
        plot(tt(startNdx_ST(i):endNdx_ST(i))./f_orco.fs,curv(startNdx_ST(i):endNdx_ST(i)),'r','LineWidth',1);
    end
    for i = 1:numel(startNdx_stop)
        plot(tt(startNdx_stop(i):endNdx_stop(i))./f_orco.fs,curv(startNdx_stop(i):endNdx_stop(i)),'g','LineWidth',1);
    end
    xlabel('time (s)');ylabel('curvature (deg/s)')
    
    % plot the firing rate of the trajectory
    subplot(4,1,4);
    plot(tt./f_orco.fs,currF,'k','LineWidth',1);hold on;
    for i = 1:numel(startNdx_ST)
        plot(tt(startNdx_ST(i):endNdx_ST(i))./f_orco.fs,currF(startNdx_ST(i):endNdx_ST(i)),'r','LineWidth',1);
    end
    for i = 1:numel(startNdx_stop)
        plot(tt(startNdx_stop(i):endNdx_stop(i))./f_orco.fs,currF(startNdx_stop(i):endNdx_stop(i)),'g','LineWidth',1);
    end
    xlabel('time (s)');ylabel('firing rate')
end

% % print the figures to pdf
% if plotFig
%     for f = 1:get(gcf,'Number')
%         figure(f);
%         print('-painters','-dpsc2',[figureFile '.ps'],'-loose','-append');
%     end
%     ps2pdf('psfile', [figureFile '.ps'], 'pdffile', ...
%         [figureFile '.pdf'], 'gspapersize', 'letter',...
%         'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
%         'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
%         'gslibpath','C:\Program Files\gs\gs9.50\lib');
% end
end