function h = plotCircle(center,r,N,color)
% plotCircle  plots a circle with a given center, radius, and color
%
% Inputs:
%    center: [x,y] location of circle center
%    r: radius of cirle
%    N: number of data points to use when plotting circle
%    color: circle color
%
% Outputs:
%    h: line property handle

theta=linspace(0,2*pi,N);
rho=ones(1,N)*r;
[X,Y] = pol2cart(theta,rho);
X=X+center(1);
Y=Y+center(2);
h=plot(X,Y,'Color',color,'LineStyle','-','Linewidth',1);