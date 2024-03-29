function [b_spdLeave,b_curvLeave,b_spdEnter,b_curvEnter] = ...
    linearFilterAnalysisSpeedCurvature(f_orco,meta,plotFigure)
% linearFilterAnalysisSpeedCurvature  wrapper function for calculating the
%   linear filters for speed and curvature at crossing (defined by delta
%   firing rate)
%
%   Inputs: f_orco = fly object
%
%   Output: b_spdLeave = n x m linear filters for speed when leaving. Each 
%               column is a filter with increasing reguarization param  
%           b_curvLeave = linear filters for curvature when leaving.
%           b_spdEnter = linear filters for speed when entering.
%           b_curvEnter = linear filters for curvature when entering.
%   
close all
if plotFigure
    if ~exist(strcat(string(meta.plotFold),'/LinearFilterAnalysis/'), 'dir')
        mkdir(strcat(string(meta.plotFold),'/LinearFilterAnalysis/'))
    end
end
fs = 100;
close all
% compute turn triggered averages
delay = 4*f_orco.fs;
buffer = 1*f_orco.fs;
%% all flies
figureFile = 'Firing Rate to Kinematics Filter Analysis 092722';
thresh = 15;
yl = 'mm/s';
[df_leaving_all,df_entering_all,f_leaving_all,f_entering_all,spd_leaving,...
    spd_entering,spd_leaving_std,spd_entering_std] = ...
    getDf_aligned_kinematics(f_orco,thresh,delay+buffer,yl,'spd',false,plotFigure);

yl = 'degrees/s';
[~,~,~,~,curv_leaving,curv_entering,curv_leaving_std,curv_entering_std] = ...
    getDf_aligned_kinematics(f_orco,thresh,delay+buffer,yl,'curv',false,plotFigure);

delay = floor(delay.*fs./f_orco.fs);
curv_leaving = resample(padarray(curv_leaving,[0 5],'replicate'),fs,f_orco.fs);
curv_leaving = curv_leaving((ceil(end./2)-delay):(ceil(end./2)+delay));
curv_entering = resample(padarray(curv_entering,[0 5],'replicate'),fs,f_orco.fs);
curv_entering = curv_entering((ceil(end./2)-delay):(ceil(end./2)+delay));
spd_leaving = resample(padarray(spd_leaving,[0 5],'replicate'),fs,f_orco.fs);
spd_leaving = spd_leaving((ceil(end./2)-delay):(ceil(end./2)+delay));
spd_entering = resample(padarray(spd_entering,[0 5],'replicate'),fs,f_orco.fs);
spd_entering = spd_entering((ceil(end./2)-delay):(ceil(end./2)+delay));

curv_leaving_std = resample(padarray(curv_leaving_std,[0 5],'replicate'),fs,f_orco.fs);
curv_leaving_std = curv_leaving_std((ceil(end./2)-delay):(ceil(end./2)+delay));
curv_entering_std = resample(padarray(curv_entering_std,[0 5],'replicate'),fs,f_orco.fs);
curv_entering_std = curv_entering_std((ceil(end./2)-delay):(ceil(end./2)+delay));
spd_leaving_std = resample(padarray(spd_leaving_std,[0 5],'replicate'),fs,f_orco.fs);
spd_leaving_std = spd_leaving_std((ceil(end./2)-delay):(ceil(end./2)+delay));
spd_entering_std = resample(padarray(spd_entering_std,[0 5],'replicate'),fs,f_orco.fs);
spd_entering_std = spd_entering_std((ceil(end./2)-delay):(ceil(end./2)+delay));

f_leaving_all = resample(padarray(f_leaving_all,[0 5],'replicate'),fs,f_orco.fs);
f_leaving_all = f_leaving_all((ceil(end./2)-delay):(ceil(end./2)+delay));
f_entering_all = resample(padarray(f_entering_all,[0 5],'replicate'),fs,f_orco.fs);
f_entering_all = f_entering_all((ceil(end./2)-delay):(ceil(end./2)+delay));


%% generate linear filters
ntfilt = 2*fs;
nthist = 0;

% firing rate filter leaving
stim_all = {f_leaving_all};filter_type = {'firing rate filter'};
[U,s,V,XStimNew] = getSVD(stim_all,ntfilt,nthist,delay);

currKin = spd_leaving(delay+1:end)';
currKin_std = spd_leaving_std(delay+1:end);
yl = 'mm/s';
b_spdLeave = generateLinearFilter_KinTurnProb(stim_all,currKin,currKin_std,ntfilt,filter_type,delay,yl,14,U,s,V,XStimNew,fs,plotFigure);
if plotFigure
    sgtitle('Speed Leaving')
end

currKin = curv_leaving(delay+1:end)';
currKin_std = curv_leaving_std(delay+1:end);
yl = 'deg/s';
b_curvLeave = generateLinearFilter_KinTurnProb(stim_all,currKin,currKin_std,ntfilt,filter_type,delay,yl,14,U,s,V,XStimNew,fs,plotFigure);
if plotFigure
    sgtitle('Curvature Leaving')
end

% firing rate filter entering
stim_all = {f_entering_all};filter_type = {'firing rate filter'};
[U,s,V,XStimNew] = getSVD(stim_all,ntfilt,nthist,delay);

currKin = spd_entering(delay+1:end)';
currKin_std = spd_entering_std(delay+1:end)';
yl = 'mm/s';
b_spdEnter = generateLinearFilter_KinTurnProb(stim_all,currKin,currKin_std,ntfilt,filter_type,delay,yl,14,U,s,V,XStimNew,fs,plotFigure);
if plotFigure
    sgtitle('Speed Entering')
end

currKin = curv_entering(delay+1:end)';
currKin_std = curv_entering_std(delay+1:end)';
yl = 'deg/s';
b_curvEnter = generateLinearFilter_KinTurnProb(stim_all,currKin,currKin_std,ntfilt,filter_type,delay,yl,14,U,s,V,XStimNew,fs,plotFigure);
if plotFigure
    sgtitle('Curvature Entering')
end

% %% individual flies
% thresh = 15;
% yl = 'mm/s';
% [df_leaving_all,df_entering_all,f_leaving_all,f_entering_all,spd_leaving,...
%     spd_entering,spd_leaving_std,spd_entering_std] = ...
%     getDf_aligned_kinematics(f_orco,thresh,delay,yl,'spd',true);
% yl = 'radians';
% [~,~,~,~,curv_leaving,curv_entering,curv_leaving_std,curv_entering_std] = ...
%     getDf_aligned_kinematics(f_orco,thresh,delay,yl,'curv',true);
% 
% %% generate linear filters
% ntfilt = 2*fs;
% nthist = 0;
% 
% % firing rate filter leaving for individual flies
% for fly = 1:f_orco.nFly
%     stim_all = f_leaving_all(fly);filter_type = {'firing rate filter'};
%     [U,s,V,XStimNew] = getSVD(f_leaving_all(fly),ntfilt,nthist,delay);
%     currKin = spd_leaving{fly}(delay+1:end)';
%     currKin_std = spd_leaving_std{fly}(delay+1:end)';
%     yl = 'mm/s';
%     b_spdLeave = generateLinearFilter_Kin(stim_all,currKin,currKin_std,ntfilt,filter_type,delay,yl,14,U,s,V,XStimNew);
%     sgtitle(['Fly ' num2str(fly) ' Speed Leaving'])
% end
% 
% %a = 1;
% % 
% % % delta firing rate filter
% % stim_all = {df_leaving_all};filter_type = {'firing rate filter'};
% % [U,s,V,XStimNew] = getSVD(stim_all,ntfilt,nthist,delay);
% % 
% % currKin = spd_leaving(delay+1:end)';
% % currKin_std = spd_leaving_std(delay+1:end);
% % generateLinearFilter_Kin(stim_all,currKin,currKin_std,ntfilt,filter_type,delay,yl,[],U,s,V,XStimNew)
% % sgtitle('Speed Leaving')
% % 
% % currKin = curv_leaving(delay+1:end)';
% % currKin_std = curv_leaving_std(delay+1:end);
% % generateLinearFilter_Kin(stim_all,currKin,currKin_std,ntfilt,filter_type,delay,yl,[],U,s,V,XStimNew)
% % sgtitle('Curvature Leaving')
% % 
% % % both filters
% % stim_all = {f_leaving_all,df_leaving_all};
% % filter_type = {'firing rate filter','delta firing rate filter'};
% % [U,s,V,XStimNew] = getSVD(stim_all,ntfilt,nthist,delay);
% % 
% % currKin = spd_leaving(delay+1:end)';
% % currKin_std = spd_leaving_std(delay+1:end);
% % generateLinearFilter_Kin(stim_all,currKin,currKin_std,ntfilt,filter_type,delay,yl,[],U,s,V,XStimNew)
% % sgtitle('Speed Leaving')
% % 
% % currKin = curv_leaving(delay+1:end)';
% % currKin_std = curv_leaving_std(delay+1:end);
% % generateLinearFilter_Kin(stim_all,currKin,currKin_std,ntfilt,filter_type,delay,yl,[],U,s,V,XStimNew)
% % sgtitle('Curvature Leaving')
% % 
% 
if plotFigure
    psFileName = strcat(string(meta.plotFold),'/LinearFilterAnalysis/',figureFile,'.ps');
    if exist(psFileName, 'file')==2
      delete(psFileName);
    end
    for f = 1:get(gcf,'Number')
        figure(f);
        print('-painters','-dpsc2',psFileName,'-loose','-append');
    end
    try
        ps2pdf('psfile', psFileName, 'pdffile', ...
        strcat(string(meta.plotFold),'/LinearFilterAnalysis/', figureFile,'.pdf'), 'gspapersize', 'letter',...
        'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
        'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
        'gslibpath','C:\Program Files\gs\gs9.50\lib');
    catch
        disp('No ghostscript available. Please install ghostscript or ')
        disp('change path to ghostscript in line 174 of linearFilterAnalysisSpeedCurvature.m')
    end
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