function [modeled_Intensity,x_tot] = Intensity_to_Space(foldernameFig,rad,Intensity)
% Intensity_to_Space takes in the radius of the fly arena and a measured
% Intensity train to use in fitting to a Intensity-distance model. The
% model is composed on a lower pass filter component and a third order
% polynomial

% Inputs:   foldernameFig: folder to save figures to
%           rad: radius of arena
%           Intensity: Measured Intensity for model fitting
% Outputs:  modeled_Intensity: New modeled intensity values
%           x_tot = Distance from edge of arena

% Note:
%   Intensity is measured and converted to nW. We want a final intensity of
%   mW/cm^2
%   You can specify delX to get finer resolution. delX is in mm
%   Plotting is currently done only for butterworth filter and not for chebyshev

%
% 2017, Liangyu Tao

% Fitting data
cF=(10^(-6))/(pi*(0.5*10^(-1))^2);                                          % conversion factor for nW to mW/cm^2 assuming radius = 0.5 mm
Intensity = Intensity*cF;
x = [0:1:length(Intensity)-1]*0.1;

%Preprocessing
ratio = Intensity(2:end)./Intensity(1:end-1);
cutoffNdx = find(ratio<0.9)-1;
Intensity_Inside = Intensity(1:cutoffNdx(1)+1);
type = 2;                                                                   %signifies using butterworth as filter

% Calculate transfer functions
transFun = Outside_TF(Intensity,type);
TF_In = Inside_TF(foldernameFig,Intensity_Inside);

%new resolution we want
delX = 1;

%specify distances
x_tot = [0:delX:rad].*0.1;
x_in = x_tot(1:((cutoffNdx(1)+1)/delX));
x_out = x_tot(((cutoffNdx(1)+1)/delX)+1:end);

% array of 1 and zero to signify if we are at the odor cutoff distance
x_bord = [ones(1,length(x_in)),zeros(1,length(x_out))];

% plot the predictions against the original measurements
if type == 2
    Wp = transFun.Wp;
    W = x_tot;
    eps = transFun.eps;
    n = transFun.n;
    
    f = (eps.^2).*(W./Wp).^(2*n);
    I = 1./sqrt(1+f);
    
    y1 = polyval(TF_In.coeff,x_tot);
    modeled_Intensity = ((y1).*Intensity(1)-max(Intensity)).*x_bord+I.*max(Intensity);
    figure
    plot(x_tot,modeled_Intensity,'b')
    hold on
    plot(x,Intensity,'r')
    hold off
    xlabel('Distance (cm)')
    ylabel('mW/cm^2')
    ylim([0 max(Intensity)*1.1])
    legend({'Modeled','Measured'})
    print('-dpdf',[foldernameFig '\Intensity_To_Distance_TF' 'Tracks.pdf']);
end


end



