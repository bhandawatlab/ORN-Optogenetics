function [transFun] = Outside_TF(Intensity,filtertype)
% Outside_TF takes in a filtertype (default is butterworth) and intensity
% measurements and creates the low pass transfer function that best fits
% the boundary and outside regions of the fly arena

% Inputs:   Intensity: Measured Intensity for model fitting
%           filtertype: type of filter. 1 = chebyshev, 2 = butterworth
% Outputs:  transFun: Low pass transfer function parameters (structure)

%
% 2017, Liangyu Tao

ratio = Intensity(2:end)./Intensity(1:end-1);
cutoffNdx = find(ratio<0.9)-1;

baseline = Intensity(end);
Intensity = Intensity - baseline;
Inside_Intensity = min(Intensity(1:cutoffNdx(1)-1));

if filtertype == 1
    %fit to chebyshev filter
    err = zeros(20,50);
    for i = 1:20
        for j = 1:50
            Wp = cutoffNdx(1)./10;
            W = 0:0.1:1.5;
            eps = 0.1+j.*0.01;
            n = i;
            
            T = cos(n*acos(W./Wp));
            f = (eps.^2).*(T.^2);
            I = 1./(1+f);
            
            modeled_Intensity = (Inside_Intensity).*I;
            
            err(i,j) = sum((modeled_Intensity(cutoffNdx(1):end)-Intensity(cutoffNdx(1):end)).^4);
        end
    end
    [bestFit(1),bestFit(2)] = find(err == min(min(err)));
    transFun.Wp = Wp;
    transFun.W = W;
    transFun.eps = 0.1+bestFit(2).*0.1;
    transFun.n = bestFit(1);
    transFun.err = min(min(err));
    transFun.type = 'Chebyshev Type 1';
    
    T = cos(transFun.n*acos(transFun.W./transFun.Wp));
    f = (transFun.eps.^2).*(T.^2);
    I = 1./(1+f);
    modeled_Intensity = (Inside_Intensity).*I;
    
    figure
    subplot(1,2,1)
    plot(W,Intensity,'r')
    hold on
    plot(W,modeled_Intensity,'b')
    hold off
    ylim([0 max(Intensity+baseline)])
    
    subplot(1,2,2)
    plot(W,Intensity+baseline,'r')
    hold on
    plot(W,modeled_Intensity+baseline,'b')
    hold off
    ylim([0 max(Intensity+baseline)])
    
end

if filtertype == 2
    %use a butterworth low pass design
    err = zeros(50,20);
    for i = 1:50
        for j = 1:20
            Wp = (cutoffNdx(1)+1)./10;
            W = 0:0.1:1.5;
            eps = 0.9+0.01*j;
            n = i;
            f = (eps.^2).*(W./Wp).^(2*n);
            I = 1./sqrt(1+f);
            
            modeled_Intensity = (Inside_Intensity).*I;
            
            err(i,j) = sum((modeled_Intensity(cutoffNdx(1):end)-Intensity(cutoffNdx(1):end)).^4);
        end
    end
    [bestFit(1),bestFit(2)] = find(err == min(min(err)));
    transFun.Wp = Wp;
    transFun.W = W;
    transFun.eps = 0.9+0.01*bestFit(2);
    transFun.n = bestFit(1);
    transFun.err = min(min(err));
    transFun.type = 'Butterworth Low Pass';
    
    Wp = (cutoffNdx(1)+1)./10;
    W = 0:0.1:1.5;
    eps = 0.9+0.01*bestFit(2);
    n = bestFit(1);
    f = (eps.^2).*(W./Wp).^(2*n);
    I = 1./sqrt(1+f);
    modeled_Intensity = (Inside_Intensity).*I;
    
    figure
    subplot(1,2,1)
    plot(W,Intensity,'r')
    hold on
    plot(W,modeled_Intensity,'b')
    hold off
    ylim([0 max(Intensity+baseline)*1.1])
    xlabel('Distance (cm)')
    ylabel('mW/cm^2')
    title('Butterworth Low Pass Filter')
    
    subplot(1,2,2)
    plot(W,Intensity+baseline,'r')
    hold on
    plot(W,modeled_Intensity+baseline,'b')
    hold off
    ylim([0 max(Intensity+baseline)*1.1])
    xlabel('Distance (cm)')
    ylabel('mW/cm^2')
end

end