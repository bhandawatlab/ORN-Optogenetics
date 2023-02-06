function [I] = Intensity_to_Time(Intensity_space,x,rad,s,foldernameFig,r)
% Intensity_to_Time converts an intensity to space model generated from the
% Intensity_to_Space.m file and a position data file to generate a
% intensity track for the position data.

% Inputs:   Intensity_space: Modeled Intensity as a function of space
%           x: distance away from the center of arena
%           rad: radius of arena in mm
%           Intensity: Measured Intensity for model fitting
%           filename = folder where the position data files are stored.
%           foldernameFig = folder where the figures are stored.
% Outputs:  I: the modeled Intensity track for the given fly position data

%
% 2017, Liangyu Tao

% Use this commented out code if we are using the old .dat files on dropbox
% % Find the radius of the arena from the circleImage of arena
% center(1) = (max(data.vertices(:,1))+min(data.vertices(:,1)))/2;
% center(2) = (max(data.vertices(:,2))+min(data.vertices(:,2)))/2;
% radius = ((max(data.vertices(:,1))-min(data.vertices(:,1)))+...
%     (max(data.vertices(:,2))-min(data.vertices(:,2))))/4;
% 
% 
% % convert position to a distance away from center of arena
% xPos = data.smoothx;
% yPos = data.smoothy;
% r = sqrt((xPos-center(1)).^2+(yPos-center(2)).^2);
% r = r./radius;

x = x./(rad/10);                                                            %x is now normalized

% Use this if we are using the new .dat files with slip, thurst, and yaw
if isempty(r)
    r = s.Distances.DistanceR;
    lightOn = s.LightOn;
else
    lightOn = 5400;
end

bin = 1:1:length(x);
I = zeros(size(r));

% convert Intensity = f(space) to f(time)
for i = 1:length(bin)-1
    ndx = r<x(i+1) & r>=x(i);
    I(ndx) = Intensity_space(i);
end
I(:,1:lightOn-1) = 0;
I(I<0.001) = 0;
%I(end-10:end) = 0;

% % Plot the Intensity trajectory for the current track
% figure
% subplot(4,1,1)
% plot(1-r,'r')
% hold on
% plot([1,10800],ones(2)*(1.7/3.2),'g')
% plot([1,10800],ones(2)*(2.2/3.2),'b')
% hold off
% xlim([0 10800])
% ylabel('1-distance')
% title('Distance from center')
% 
% subplot(4,1,2)
% plot(1-r,'r')
% hold on
% plot([1,10800],ones(2)*(2.2/3.2),'b')
% hold off
% ylim([1.7/3.2,1])
% xlim([0 10800])
% ylabel('1-distance')
% title('Distance from center cp = 1.5')
% 
% subplot(4,1,3)
% plot(1-r,'r')
% ylim([2.2/3.2,1])
% xlim([0 10800])
% ylabel('1-distance')
% title('Distance from center cp = 1')
% 
% subplot(4,1,4)
% ylabel('1-distance')
% plot(I,'r')
% xlim([0 10800])
% ylabel('mW/cm^2')
% set(gcf,'position',[9    49   944   948])
%print('-dpdf',[foldernameFig '\' filename(50:end-8) 'Tracks.pdf']);

end