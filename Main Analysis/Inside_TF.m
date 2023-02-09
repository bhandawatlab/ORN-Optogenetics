function [transFun_Inside] = Inside_TF(foldernameFig,Intensity)
% Inside_TF takes in the inside arena intensity values and fits a 3rd order
% polynomial to it.

% Inputs:   foldernameFig: folder to save figures to
%           Intensity: Measured Intensity for model fitting
% Outputs:  transFun_Inside: 3rd order polynomial parameters (structure)

% Note:
%   Intensity is only of the inside stim ring intensities
%   A higher order polynomial can by specified by changing the third input
%   to polyit

%
% 2017, Liangyu Tao

ratio = Intensity/Intensity(1);
x = [1:1:length(Intensity)]*0.1;
ratio = smooth(ratio);

p = polyfit(x,ratio',3);
y1 = polyval(p,x);

modeled_Intensity = y1.*Intensity(1);

plot(x,Intensity,'r');
hold on
plot(x,modeled_Intensity,'b')
ylim([0 max(Intensity)*1.1])
xlabel('Distance (cm)')
ylabel('mW/cm^2')
title('3rd Order Polynomial')
print('-dpdf',[foldernameFig '\I_to_D_TFComposition' 'Tracks.pdf']);

transFun_Inside.coeff = p;

end