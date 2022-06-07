function [] = plotXYTracks(f_orco,border,nFly,sort_by_attraction,id)

fe = f_orco.getFirstEntry('H',border);
attn = f_orco.getAttnNdx('H',border);
attn = attn.during;
if sort_by_attraction
    [~,ndx] = sort(attn);
    x = f_orco.x(ndx,:);
    y = f_orco.y(ndx,:);
    fe = fe(ndx);
    attn = attn(ndx);
else
    x = f_orco.x;
    y = f_orco.y;
end

if isempty(nFly)
    nFly = f_orco.nFly;
end

for fly = 1:min(nFly,f_orco.nFly)
    if mod(fly,20)==1
        figure;set(gcf,'Position',[2 42 838 924])
        suptitle([f_orco.id ' ' id])
        k = 1;
    end
    
    subplot(5,4,k);hold on;
    plot(x(fly,1:fe(fly)),y(fly,1:fe(fly)),'g');hold on;
    plot(x(fly,fe(fly):end),y(fly,fe(fly):end),'r');
    plotCircle([0 0],4,100,'k');
    plotCircle([0 0],border,100,'b');
    xlim([-4 4]);ylim([-4 4])
    axis square
    title(['Attr Ndx=' num2str(round(attn(fly),3))])
    k = k+1;
end
    
end