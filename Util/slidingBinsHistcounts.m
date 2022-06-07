function [N,xCent,yCent,xRange,yRange,raw] = slidingBinsHistcounts...
    (x,y,dx,dy,xrange,yrange,nOvlapX,nOvlapY,histcounts)

xx = xrange(1):(dx-nOvlapX):xrange(2);
yy = yrange(1):(dy-nOvlapY):yrange(2);


if histcounts == false
    xRange = [xx(1:end-1);xx(1:end-1)+dx];
    yRange = [yy(1:end-1);yy(1:end-1)+dy];
    N = nan(numel(xx)-1,numel(yy)-1);
    raw = cell(numel(xx)-1,numel(yy)-1);
else
    N = zeros(numel(xx)-1,numel(yy)-1);
    xRange = zeros(2, numel(xx)-1);
    yRange = zeros(2, numel(yy)-1);
    raw = cell(numel(xx)-1,numel(yy)-1);
    for i = 1:numel(xx)-1
        for j = 1:numel(yy)-1
            mask = (x>xx(i) & x<=(xx(i)+dx)) & (y>yy(j) & y<=(yy(j)+dy));
            raw{i,j} = [x(mask), y(mask)];
            N(i,j) = sum(mask(:));
            xRange(:,i) = [xx(i); xx(i)+dx];
            yRange(:,j) = [yy(j); yy(j)+dy];
        end
    end

end

xCent = xx(1:end-1)+dx/2;
yCent = yy(1:end-1)+dy/2;

end