function [synth_orco,f_orco] = GenerateSyntheticFlies(fName,gen,meta,savefile)
C = strsplit(fName,'.');
fName2 = [C{1} '_flies.mat'];
[synth_orco,f_orco] = getSyntheticFlies([meta.syntheticFlyFold fName],gen,meta);
if savefile
    save([meta.syntheticFlyFold fName2],'synth_orco','f_orco');
end

end

function [synth_orco,f_orco] = getSyntheticFlies(fName,gen,meta)

%load(['DataRT/RT_run' gen '_' meta.d cond '.mat'],'-v7.3');
load(fName,'synthFlys','params','state','curvAll','dfSmooth','f_orco','LLFfs');
fs = LLFfs;
rBound = meta.rBound;
border = meta.border;

% get the baseline
baseline = mean(synthFlys.spk(:,30*fs:60*fs),'all');
% remove any artifacts from the linear filter
synthFlys.spk(:,1:30*fs) = round(baseline,4);
synthFlys.spk = round(synthFlys.spk,4);
synthFlys.spk = synthFlys.spk(:,1:size(synthFlys.x,2));

% create flies object
synth_orco = Flies(gen,synthFlys.x,synthFlys.y,synthFlys.x,synthFlys.y,...
    synthFlys.spk,fs,rBound,synthFlys.lightOn,[],[],[],[],size(synthFlys.x,2));

synth_orco.states.ndx = state;
synth_orco.states.key = {'sharp turns';'curved walks';'stops';'boundary'};
synth_orco.curv = curvAll;

% calculate first entry
fe = synth_orco.getFirstEntry('H',border);

% remove flies that do not enter
synth_orco = synth_orco.rmvData(isnan(fe));
fe = synth_orco.getFirstEntry('H',border);

% get the firing rate and change in firing rate
spk = synth_orco.spk;
dSpk = synth_orco.calcDeltaFR;

% separate out into each time point into before first entry (FE), below
% baseline firing rate, baseline firing rate after FE, and above
% baseline firing rate
condNdx = zeros(size(spk));condNdx(spk>baseline) = 4;
condNdx(spk<baseline) = 2;condNdx(abs((spk-baseline))<0.001 & dSpk==0) = 3;
for j = 1:synth_orco.nFly
    condNdx(j,1:fe(j)-1) = 1;
end
key = {'before','below','baseline','above'};

meta.zGrid = meta.zGrid.*synth_orco.fs./f_orco.fs;
%meta.ratio = f_orco.model.params{1,1}.KNN.ratio;
meta.ratio(3) = meta.ratio(3).*synth_orco.fs./f_orco.fs;

% get the kinematics and decision space of the flies
GetKinematicModelParams(synth_orco,condNdx,key,fe,meta,true);

% intensity space of the arena
load('Data/Intensity_space2.mat','xN','p','Intensity_spaceN','convIV','convVI')
synth_orco.model.IntensitySpace.x = xN;
synth_orco.model.IntensitySpace.I = Intensity_spaceN;
synth_orco.model.IntensitySpace.p = p;
synth_orco.model.IntensitySpace.convIV = convIV;% conversion from intensity to voltage
synth_orco.model.IntensitySpace.convIV = convVI;% conversion from voltage to intensity

% when df is 0 and there is a nonbaseline firing rate (inhibition period)
%crossingtype = 'enter';
crossingtype = 'exit';%plotFig = false;
[synth_orco,~,~,~,~,~,~,~] = getDistInhibition(synth_orco,crossingtype,meta.plotFig);
end





















