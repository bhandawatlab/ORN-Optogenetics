function [I,V] = GenerateLightIntensity(label,meta)
% light_Calibration_Wrapper is the main wrapper function for converting a
% track of the flies centroid into voltage traces for manipulating the
% light intensity in optogenetics electrophys.

% Inputs: label = '' or 'Head'
%
% Outputs: I = Intensity trace in mW/cm^2
%          V = Voltage trace matching to the intensity

fs = meta.fs_linFilter;
fs_beh = meta.fs;
foldernameFig = [meta.DestPath '\Figures'];
rad = meta.rBound*10;% in mm
Intensity = [26030, 26020, 26090, 25630, 25800, 25440, 25500, 24920, ...    % Measured intensity in nW
    25030, 24780, 25100, 16700, 3610, 134, 31.5, 20];

% Transfer Function model Fitting
[Intensity_space,x] = Intensity_to_Space(foldernameFig,rad,Intensity);
% upsample by 100x, so that intensity jumps aren't as jagged
xN = x(1):mean(diff(x))/100:x(end);
Intensity_spaceN = interp1(x,Intensity_space,xN);

dirinfo = dir(string(meta.folderData));
subdirinfo = cell(1,length(dirinfo)-2);
for K = 1 : length(subdirinfo)
    thisdir = dirinfo(K+2).name
    %[string(meta.folderData) '/' thisdir]
    subdirinfo{K} = dir(strcat(string(meta.folderData),'/',thisdir));
end

% Converting position data to Voltage traces
[p,~] = Voltage_to_Intensity_Converter(string(meta.calibrationFolder),string(meta.calibrationFile));
convIV = @(I) (round((I-p(2))/p(1),3)>p(3)).*round((I-p(2))/p(1),3);
convVI = @(V) max(round(p(1).*V+p(2),3),0);
save(string(meta.Intensity2VoltageFile),'Intensity_spaceN','xN','convIV','convVI','p','-v7.3')

for i = 1:length(subdirinfo)
    close all
    filename = strcat(string(meta.folderData),'\',subdirinfo{i}.name);
    load(filename,'Data');
    r = sqrt(Data.(['x' label]).^2 + Data.(['y' label]).^2)./4;r(r>1) = 1;
    %r = empFlys.rH./4;
    [I] = Intensity_to_Time(Intensity_spaceN,xN,rad,filename,foldernameFig,r);
    V = (I-p(2))/p(1);
    V = round(V,3);
    p(3) = min(V(:));
    V(V==p(3)) = 0;
    V2 = V;
    I2 = I;
    C = strsplit(subdirinfo{i}.name,'_');
    if ~isempty(V2)
        % Decrease the sampling rate from 30 Hz to 100 Hz
        V2 = resample(V2',fs,fs_beh,0)';
        I2 = resample(I2',fs,fs_beh,0)';
        matfile = strcat(string(meta.foldStim),'\',C{1},'_Stim_Train.mat');
        save(matfile, '-v7.3','I','I2','V','V2','convIV','convVI','p','fs');
    end
    if isempty(V2)
        matfile = strcat(string(meta.foldStim),'\Bad_Data\',C{1},'_Stim_Train.mat');
        save(matfile, '-v7.3','I','I2','V','V2','convIV','convVI','p','fs');
    end
end



end