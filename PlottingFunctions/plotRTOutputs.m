function [] = plotRTOutputs(synth_orco,f_orco,border)
fNum = get(gcf,'Number');

% get some basic stuff like attraction index, radial probability, and
% decision density, and probability of being inside
fe = synth_orco.getFirstEntry('H',border);
attNdx = synth_orco.getAttnNdx('H',border);
radProb = synth_orco.getRadialProb('H',border,0,synth_orco.fs);
%turnDens = synth_orco.getTurnDens('H',s,keyNdx);
probIn = synth_orco.getProbIn('H',border,fe);
probIn2 = synth_orco.getProbIn('H',border,ones(synth_orco.nFly,1));

feEmp = f_orco.getFirstEntry('H',border);
attNdxEmp = f_orco.getAttnNdx('H',border);
radProbEmp = f_orco.getRadialProb('H',border,0,f_orco.fs);
%turnDens = f_orco.getTurnDens('H',s,keyNdx);
probInEmp = f_orco.getProbIn('H',border,feEmp);
probIn2Emp = f_orco.getProbIn('H',border,ones(f_orco.nFly,1));

% badNdx = attNdx.duringTot<1000;
% attNdx.during(badNdx) = nan;
[sortedAttnNdx,ndx] = sort(attNdx.during);
[~,midNdx] = min(abs(sortedAttnNdx-0.5));%floor((synth_orco.nFly-sum(badNdx))./2);

[~,ndx2] = sort(attNdxEmp.during);

%flys2Cons = [ndx(1:10) ndx(midNdx-5:midNdx+4) ndx(end-sum(badNdx)-9:end-sum(badNdx))];
flys2Cons = [ndx(1:15) ndx(end-14:end)];

figure(fNum+3);set(gcf,'Position',[1 40 585 600])
for i = 1:numel(flys2Cons)
    fNdx = flys2Cons(i);
    subplot(6,5,i);
    plotCircle([0 0],border,100,'k');hold on
    plot(synth_orco.x(fNdx,1:fe(fNdx)),synth_orco.y(fNdx,1:fe(fNdx)),'g');
    plot(synth_orco.x(fNdx,fe(fNdx):end),synth_orco.y(fNdx,fe(fNdx):end),'r');hold off
    axis([-4 4 -4 4])
    title(['Attn Ndx = ' num2str(attNdx.during(fNdx))])
end

p = tempPlottingFunc(f_orco,attNdxEmp,radProbEmp,probInEmp,probIn2Emp,border,1,[],fNum+1);
tempPlottingFunc(synth_orco,attNdx,radProb,probIn,probIn2,border,2,p,fNum+1);



subplots = [2,1,1,fNum+2]; % ysbplt, xsbplt, sbpltN, fig
xSubplots = repmat(100./subplots(1),1,subplots(1));
ySubplots = repmat(100./subplots(2),1,subplots(2));
figure(subplots(4));set(gcf,'Position',[1 40 585 600])
p = panel();
p.pack(xSubplots, ySubplots);
subplots(3) = 1;%1

% define basic paramters
lims = [0 1];                      % y limits to set (example: [0 1])
isPaired = 'N';                 % equivalent to paired t-test
circleSize = 60;                % size of scatter plot points
barstate = 'off';                % 'On' = bar graph, 'Off' = scatter plot

dat1 = [attNdxEmp.before attNdx.before];
g1 = [repmat({'emp'},f_orco.nFly,1); repmat({'synth'},synth_orco.nFly,1)];
[ss,p,~] = dabest3(dat1,g1,p,[],[0 1],isPaired,circleSize,barstate,subplots);
% set axis labels for comparisons against category 1
[J,I] = ind2sub(subplots(2:-1:1),subplots(3));
try
    p(I,J).select();
    ylabel('atn Ndx');
    title([synth_orco.id ' before'],'Interpreter','none')
catch
    p(I,J,1,1).select();
    ylabel('atn Ndx');
    title(synth_orco.id,'Interpreter','none')
    p(I,J,2,1).select();
    ylabel('delta atn Ndx');
    xtickangle(10)
end
subplots(3) = subplots(3)+1;

dat1 = [attNdxEmp.during attNdx.during];
g1 = [repmat({'emp'},f_orco.nFly,1); repmat({'synth'},synth_orco.nFly,1)];
[ss,p,~] = dabest3(dat1,g1,p,[],[0 1],isPaired,circleSize,barstate,subplots);
% set axis labels for comparisons against category 1
[J,I] = ind2sub(subplots(2:-1:1),subplots(3));
try
    p(I,J).select();
    ylabel('atn Ndx');
    title([synth_orco.id ' during'],'Interpreter','none')
catch 
    p(I,J,1,1).select();
    ylabel('atn Ndx');
    title(synth_orco.id,'Interpreter','none')
    p(I,J,2,1).select();
    ylabel('delta atn Ndx');
    xtickangle(10)
end
subplots(3) = subplots(3)+1;

end