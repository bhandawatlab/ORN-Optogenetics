function [] = GenerateSpikingData(meta)
addpath(genpath([pwd '/Subroutines']))

% load([meta.foldStim meta.LFPFilterFile],'b','optFitN');bLFP = b(:,optFitN);
% load([meta.foldStim meta.RateFilterFile],'b','optFitN');bSpkRate = b(:,optFitN);
load(meta.LFPFilterFile,'b','optFitN');bLFP = b(:,optFitN);
load(meta.RateFilterFile,'b','optFitN');bSpkRate = b(:,optFitN);
fs_filter = 100;R = [];
fs = 30;

progressbar(0,0)
s = dir([meta.foldStim '\*.mat']);tic
for i = 1:numel(s)
    C = strsplit(s(i).name,'_');
    gen = C{1};
    
    load([meta.foldStim '\' s(i).name],'V2','V')
    ncells = size(V2,1);nPts = size(V2,2);
    nTrials_part = 2;
    nPartitions = ceil(ncells./nTrials_part);
    LFP_fit = zeros(nPts,ncells);
    sps_pred = zeros(ncells,nPts);
    sps_pred2 = zeros(ncells,ceil(nPts*30/fs_filter));
    for p = 1:nPartitions
        currP = (p-1)*nTrials_part+1:min(p*nTrials_part,ncells);
        v2_tmp = V2(currP,:);
        [LFP_fit(:,currP),ntFilt] = GenLFP(v2_tmp,bLFP,[],size(v2_tmp,1),fs_filter);%toc
        sps_pred(currP,:) = Lin_Filter(LFP_fit(:,currP),R,bSpkRate,size(v2_tmp,1),ntFilt,fs_filter);%toc
        sps_pred2(currP,:) = resample(sps_pred(currP,:)',fs,fs_filter,0)';%toc
        progressbar([],p/nPartitions)
    end
    save([meta.foldSpk '\' gen '_SpkRate.mat'],'sps_pred2','sps_pred','LFP_fit','ntFilt','V')
    progressbar(i/numel(s))
end

end