function [] = GenerateSpikingData(genAll,meta)
% GenerateSpikingData  Computes the firing rate by convolving a 2-stage 
%   linear filter with the input light intensity
%
%   Inputs: meta = a structure with fields labeling where the stimulation
%   and linear filter files are. Also includes a field of where to save the
%   firing rate data
%   

addpath(genpath([pwd '/Subroutines']))

% load in the LFP and firing rate filters
load(string(meta.LFPFilterFile),'b','optFitN');bLFP = b(:,optFitN);
load(string(meta.RateFilterFile),'b','optFitN');bSpkRate = b(:,optFitN);
fs_filter = 100;% filter sampling rate (Hz)
fs = 30;% video sampling rate (Hz)

progressbar(0,0)
for i = 1:numel(genAll)
    gen = genAll{i};
    
    % use the head position files and ignore the body position files
    %if ~contains(s(i).name,'body','IgnoreCase',true)
        load(strcat(string(meta.foldStim),'\',gen,'_Stim_Train.mat'),'V2','V')
        maxV = max(V2(:));
        insideSS = calculateSS(maxV,bLFP,bSpkRate,fs_filter);
        baseline = calculateSS(0,bLFP,bSpkRate,fs_filter);
        
        % set up partitions and preallocate matrices
        ncells = size(V2,1);nPts = size(V2,2);
        nTrials_part = 2;
        nPartitions = ceil(ncells./nTrials_part);
        LFP_fit = zeros(nPts,ncells);
        sps_pred = zeros(ncells,nPts);
        sps_pred2 = zeros(ncells,ceil(nPts*30/fs_filter));
        for p = 1:nPartitions
            currP = (p-1)*nTrials_part+1:min(p*nTrials_part,ncells);
            v2_tmp = V2(currP,:);
            % calculate the LFP and then the firing rate
            [LFP_fit(:,currP),ntFilt] = GenLFP(v2_tmp,bLFP,[],size(v2_tmp,1),fs_filter);%toc
            sps_pred(currP,:) = genFiringRate(LFP_fit(:,currP),bSpkRate,size(v2_tmp,1),ntFilt);%toc
            % downsampe the firing rate to the video sampling rate
            sps_pred2(currP,:) = resample(sps_pred(currP,:)',fs,fs_filter,0)';%toc
            progressbar([],p/nPartitions)
        end
        % save the firing rate data
        save(strcat(string(meta.foldSpk),'\',gen,'_SpkRate.mat'),'sps_pred2','sps_pred','LFP_fit','ntFilt','V','insideSS','baseline')
    %end
    progressbar(i/numel(genAll))
end

end


function [insideSS] = calculateSS(maxV,bLFP,bSpkRate,fs_filter)
totalFiltDur = size(bLFP,1)+size(bSpkRate,1);
samp = maxV.*[zeros(1,fs_filter.*2),ones(1,totalFiltDur.*2),zeros(1,fs_filter.*2)];
[LFP_fit,ntFilt] = GenLFP(samp,bLFP,[],size(samp,1),fs_filter);%toc
sps_pred = genFiringRate(LFP_fit,bSpkRate,size(samp,1),ntFilt);%toc
if maxV>0
    sps_pred = sps_pred(samp>0);
end
insideSS = mean(sps_pred(gradient(sps_pred)==0));
end

