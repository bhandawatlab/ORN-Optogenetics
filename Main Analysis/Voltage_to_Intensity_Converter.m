function [p,Rsq] = Voltage_to_Intensity_Converter(folder,file)
% Voltage_to_Intensity_Converter calibrates the relationship between volage
% input to a LED and the Intensity output. 

% Inputs:   folder: Folder where the calibration files are located
%           file: Voltage_Intensity calibration file
% Outputs:  p: linear model fit
%           Rsq: R-squared value for linear fit

% Note:
%   delta is specified when making the voltage intensity calibration file
%   All wattages are first converted to nW before analysis is conducted

%
% 2017, Liangyu Tao

% Processing of calibration file
delta = 1/3;
disp([folder file])
[~,~,raw] = xlsread([strcat(folder,file)]);
wattage = raw(:,2);
wattage = wattage(2:end);
wattageVal = cellfun(@(x) x(1:end-1), wattage, 'un', 0);
wattageUnit = cellfun(@(x) x(end), wattage, 'un', 0);
ndx = strfind(wattageUnit, 'u');
ndxU = find(not(cellfun('isempty', ndx)));
wattageVal = str2num(cell2mat(wattageVal));
wattageVal(ndxU) = wattageVal(ndxU)*(10^3);

ndxChange = [find(abs(diff(wattageVal))>100)];
ndxChange(diff(ndxChange)<5) = [];
AvgIntensity = zeros(10,1);

% Convert from Intensity to a area Intensity
cF=(10^(-6))/(pi*(0.5*10^(-1))^2);                                          % conversion factor for nW to mW/cm^2 assuming radius = 0.5 mm
wattageDensity = wattageVal*cF;

% Take average intensity
for i = 1:length(ndxChange)-1
    AvgIntensity(i,1) = mean(wattageDensity(ndxChange(i)+1:ndxChange(i+1)-1));
end

% Linear regression fit and compute R-squared
x = delta:delta:delta*10;
[p,~] = polyfit(x,AvgIntensity',1);
model = p(1)*x+p(2);
Rsq = 1 - sum((AvgIntensity' - model).^2)/sum((AvgIntensity' - mean(AvgIntensity')).^2);

% Plot fit of linear regression
figure
subplot(211)
plot(wattageDensity)
xlabel('time');
ylabel('mW/cm^2');
subplot(212)
plot(x,AvgIntensity','g*')
hold on
plot(x,model','r-')
hold off
xlabel('voltage');
ylabel('mW/cm^2');
legend('Measured','Linear Model')
set(gcf,'position',[9    49   944   948])

end