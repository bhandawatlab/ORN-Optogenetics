function [] = agentModelHandle(genAll,meta,agentModelMeta,generateFlies)

delay = agentModelMeta.delay;
dur = agentModelMeta.dur;


for g = 1:numel(genAll)
    gen = genAll{g};
    for ii = 1:numel(dur)
        for jj = 1:numel(delay)
            currDur = dur(ii);
            currDelay = delay(jj);
            cond = ['_CWdur' strrep(num2str(currDur),'.','') ...
                '_CWdelay' strrep(num2str(currDelay),'.','')];
            
            if generateFlies
                progressbar(0,0,0)
                load(strcat(string(meta.foldDataModel),'\',gen,'_',meta.d,meta.ext,'.mat'),'f_orco');
                % load information for calculating spike rate
                LLFfs = 100; % linear filters for Intensity to ORN was calculated at 100 Hz
                load(string(meta.LFPFilterFile),'optFitN','b');
                bLFP = b(:,optFitN);
                load(string(meta.RateFilterFile),'optFitN','b');
                bRate = b(:,optFitN);
                
                disp("Here 1")
                % set up synthetic fly model parameters
                params.fs = f_orco.fs;%30 hz
                params.len = f_orco.nPt;%1080;
                params.totFly = 15;
                params.posInit = [];
                params.border = meta.border./meta.rBound;
                params.lightOnTime = ceil(params.len./2);%540;% time of light on
                params.bLFP = bLFP;
                params.bRate = bRate;
                params.LLFfs = LLFfs;
                params.delay = currDelay;
                params.dur = currDur;
                nIt = 10;
                params.baseLineFR = f_orco.spk(1);
                
                %RunAndTumbleFinal032722(f_orco,params);
                %RunAndTumbleKNNLinFilter(f_orco,params);
                disp("Here 2")
                p = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(p)
                    p = parpool('local');
                end
                tic;
                parfor i = 1:nIt
                    %i = 1;
                    [synthFlysC{i},stateC{i},nTurnC{i},spdAllC{i},curvAllC{i},...
                        phiAllC{i},dSpkC{i},allAngC{i},lambda_CWC{i},lambda_StopC{i},...
                        dfSmoothC{i},fSmoothC{i},tt_baseline{i}]...
                        = RunAndTumbleKNNLinFilter(f_orco,params);
                end
                disp("Here 3")
                toc;
                synthFlys = synthFlysC{1};
                %------
                for i = 2:nIt
                    synthFlys.x = [synthFlys.x; synthFlysC{i}.x];
                    synthFlys.y= [synthFlys.y; synthFlysC{i}.y];
                    synthFlys.r = [synthFlys.r; synthFlysC{i}.r];
                    synthFlys.firstEntry = [synthFlys.firstEntry; synthFlysC{i}.firstEntry];
                    synthFlys.lightOn = [synthFlys.lightOn; synthFlysC{i}.lightOn];
                    synthFlys.spk = [synthFlys.spk; synthFlysC{i}.spk];
                end
                %------
                disp("Here 4")
                state = cell2mat(stateC');
                nTurn = cell2mat(nTurnC');
                spdAll = cell2mat(spdAllC');
                curvAll = cell2mat(curvAllC');
                phiAll = cell2mat(phiAllC');
                dSpk = cell2mat(dSpkC');
                allAng = cell2mat(allAngC');
                lambda_CW = cell2mat(lambda_CWC');
                lambda_Stop = cell2mat(lambda_StopC');
                dfSmooth = cell2mat(dfSmoothC');
                fSmooth = cell2mat(fSmoothC');
                %------
                
                % calculate the different states
                synthFlys2 = synthFlys;
                
                % remove any flies with nan's in their track
                badTracks = any(isnan(synthFlys.r),2);
                synthFlys.r(badTracks,:) = [];
                synthFlys.x(badTracks,:) = [];
                synthFlys.y(badTracks,:) = [];
                synthFlys.spk(badTracks,:) = [];
                synthFlys.lightOn(badTracks,:) = [];
                state(badTracks,:) = [];
                spdAll(badTracks,:) = [];
                curvAll(badTracks,:) = [];
                phiAll(badTracks,:) = [];
                lambda_CW(badTracks,:) = [];
                lambda_Stop(badTracks,:) = [];
                dfSmooth(badTracks,:) = [];
                synthFlys.firstEntry(badTracks,:) = [];
                
                % keep only synthetic flies that enter within the 95th quantile
                % after light on
                fe = f_orco.getFirstEntry('H',meta.border);
                feEmp = floor(quantile(fe,0.85));%95
                feThresh = ceil(feEmp.*params.LLFfs./params.fs);%8500;
                synthFlys.r(synthFlys.firstEntry>feThresh,:) = [];
                synthFlys.x(synthFlys.firstEntry>feThresh,:) = [];
                synthFlys.y(synthFlys.firstEntry>feThresh,:) = [];
                synthFlys.spk(synthFlys.firstEntry>feThresh,:) = [];
                synthFlys.lightOn(synthFlys.firstEntry>feThresh,:) = [];
                state(synthFlys.firstEntry>feThresh,:) = [];
                spdAll(synthFlys.firstEntry>feThresh,:) = [];
                curvAll(synthFlys.firstEntry>feThresh,:) = [];
                phiAll(synthFlys.firstEntry>feThresh,:) = [];
                lambda_CW(synthFlys.firstEntry>feThresh,:) = [];
                lambda_Stop(synthFlys.firstEntry>feThresh,:) = [];
                dfSmooth(synthFlys.firstEntry>feThresh,:) = [];
                synthFlys.firstEntry(synthFlys.firstEntry>feThresh,:) = [];
                
                fName = strcat('RT_run',gen,'_',meta.d,meta.ext,cond,'.mat');
                fName2 = strcat('RT_run',gen,'_',meta.d,meta.ext,cond,'_flies.mat');
                
                % save the data
                save(strcat(string(meta.syntheticFlyFold),string(fName)),'-v7.3');%_BCNew2
                disp("Here 5")
                [synth_orco,f_orco] = GenerateSyntheticFlies(string(fName), gen,meta,true);
                toc;
                cab(1);
                progressbar(g./numel(genAll),ii./numel(dur),jj./numel(delay))
            else
                fName = strcat('RT_run',gen,'_',meta.d,cond,'.mat');%['RT_run' gen '_' meta.d '' cond '.mat'];
                fName2 = strcat('RT_run',gen,'_',meta.d,meta.ext,cond,'_flies.mat');%['RT_run' gen '_' meta.d meta.ext cond '_flies.mat'];
                [synth_orco,f_orco] = GenerateSyntheticFlies(string(fName), gen,meta,false);
                save(strcat(string(meta.syntheticFlyFold),fName2),'synth_orco','f_orco');
            end
            %disp(strcat(string(meta.foldDataModel),'\',gen,'_',meta.d,meta.ext,'.mat'))
            %end
        end
    end
end

end