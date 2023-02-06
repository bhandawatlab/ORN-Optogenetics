function [theta] = get2VecAngle(v,u,n,plotFig)
% get2VecAngle  Find the angle between 2 vectors
%
% Inputs:
%    v: vector 1
%    u: vector 2
%    n: vector normal to plane that vector 1 and 2 resides in

% Outputs:
%    theta: angle in degrees between bectors 1 and 2


if isempty(n)
    theta = acosd(dot(v,u)./(norm(v).*norm(u)));
    u = u./norm(u);
    v = v./norm(v);
else
    theta = atan2d(dot(n,cross(v,u)),dot(v,u));
    u = u./norm(u);
    v = v./norm(v);
end


if plotFig
    figure;
    quiver(-u(1),-u(2),u(1),u(2));hold on
    quiver(-u(1),-u(2),v(1),v(2));
    axis equal
    legend({'u', 'v'});title([num2str(theta) ' deg'])
end
end