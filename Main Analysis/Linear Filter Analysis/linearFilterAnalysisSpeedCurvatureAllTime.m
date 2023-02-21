function [b_spd,b_curv] = linearFilterAnalysisSpeedCurvatureAllTime(f_orco,meta)
% linearFilterAnalysisSpeedCurvatureAllTime  wrapper function for calculating
%   the linear filters for speed and curvature across the entire after
%   first entry period for each fly separately
%
%   Inputs: f_orco = fly object
%           border = radial distance of light border (in cm)
%
%   Output: b_spd = fly x 1 cell array of linear filters from firing rate 
%               to speed. Each cell is for a different fly and is a matrix
%               where each column is a different filter sorted by 
%               regularization parameter lambda
%           b_curv = same as b_spd, but for linear filters from firing rate
%               to curvature.
%   
close all
if ~exist(strcat(string(meta.plotFold),'/LinearFilterAnalysis/'), 'dir')
    mkdir(strcat(string(meta.plotFold),'/LinearFilterAnalysis/'));
end

close all
% compute turn triggered averages
delay = 5*f_orco.fs;
ntfilt = 2*f_orco.fs;
%% all flies
figureFileCurv = 'Firing Rate to Curv Filter All Flies_2sFilter';
figureFileSpd = 'Firing Rate to Speed Filter All Flies_2sFilter';

fe = f_orco.getFirstEntry('H',meta.border);
spd = f_orco.calcSpd;
XStimAllFly = [];
curvAllfly = [];
spdAllfly = [];
fAllfly = [];
b_curv = cell(f_orco.nFly,1);
b_spd = cell(f_orco.nFly,1);
filter_type = {'single fly filter'};
curvPSFileName = strcat(string(meta.plotFold),'/LinearFilterAnalysis/', figureFileCurv,'.ps');
if exist(curvPSFileName, 'file')==2
    delete(curvPSFileName);
end
spdPSFileName = strcat(string(meta.plotFold),'/LinearFilterAnalysis/',figureFileSpd,'.ps');
if exist(spdPSFileName, 'file')==2
    delete(spdPSFileName);
end
for fly = 1:f_orco.nFly
    curv = abs(f_orco.curv(fly,fe(fly):end)).*180./pi.*f_orco.fs;
    curv = smooth(curv,'moving',0.2*f_orco.fs);
    currSpd = abs(spd(fly,fe(fly):end))';
    
    f = f_orco.spk(fly,fe(fly)-delay:end);
    [U,s,V,XStimNew] = getSVD({f},ntfilt,0,delay);
    
    
    yl = 'degrees/s';
    try
        b_curv{fly} = generateLinearFilter_KinAllTime({f},curv,...
            zeros(size(curv)),ntfilt,filter_type,yl,14,U,s,V,XStimNew,f_orco.fs);%delay
        sgtitle(['Fly ' num2str(fly)]);
    catch
        figure;title(['Not enough data for fly ' num2str(fly)]);
    end
    print('-painters','-dpsc2',curvPSFileName,'-loose','-append');
    yl = 'mm/s';
    try
        b_spd{fly} = generateLinearFilter_KinAllTime({f},currSpd,...
            zeros(size(currSpd)),ntfilt,filter_type,yl,14,U,s,V,XStimNew,f_orco.fs);%delay
        sgtitle(['Fly ' num2str(fly)]);
    catch
        figure;title(['Not enough data for fly ' num2str(fly)]);
    end
    print('-painters','-dpsc2',spdPSFileName,'-loose','-append');
    
    XStimAllFly = [XStimAllFly;XStimNew];
    curvAllfly = [curvAllfly;curv];
    spdAllfly = [spdAllfly;currSpd];
    fAllfly = [fAllfly,f];
end

try
    ps2pdf('psfile',curvPSFileName, 'pdffile', ...
        strcat(string(meta.plotFold),'/LinearFilterAnalysis/',figureFileCurv,'.pdf'), 'gspapersize', 'letter',...
        'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
        'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
        'gslibpath','C:\Program Files\gs\gs9.50\lib');
    ps2pdf('psfile',spdPSFileName, 'pdffile', ...
        strcat(string(meta.plotFold),'/LinearFilterAnalysis/',figureFileSpd,'.pdf'), 'gspapersize', 'letter',...
        'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
        'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
        'gslibpath','C:\Program Files\gs\gs9.50\lib');
catch
    disp('No ghostscript available. Please install ghostscript or ')
    disp('change path to ghostscript in line 81 of linearFilterAnalysisSpeedCurvatureAllTime.m')
end

end

function [U,s,V,XStimNew] = getSVD(stim_all,ntfilt,nthist,delay)
opts.ncells = 1;
XStimNew_all = [];
for i = 1:numel(stim_all)
    stim = stim_all{i};
    % design matrix
    [XStimNew,~,ntfilt,nthist] = createDesignMat(stim,[],ntfilt,nthist,opts);
    if i == 1
        XStimNew_all = XStimNew;
    else
        XStimNew_all = [XStimNew_all,XStimNew(:,2:end)];
    end
end
XStimNew = XStimNew_all(delay+1:end,:);

% prep for the ridge regression using svd
tic;[U,s,V] = csvd(XStimNew,[]);toc
end