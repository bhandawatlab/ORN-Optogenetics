clear all; close all; daqreset;
imaqreset;
%% Overall, we want to acquire a 6 minute video. During the last three minutes we will
%stimulate with red light to activate ORNs. We have to acquire video and in
%sync stimulate.
% basic parameters
%numSecs=30;
numSecs=360;
numafter=40;
sampleRate=1000;
frameRate=30;
numFrames=numSecs*frameRate;
extshutter=18000; %12000;
brightness=63; %110;
gain = 1;

%define where things will be saved
rootdir = ['F:\Liangyu Tao'];
mkdir(rootdir);
cd(rootdir);

%update the expNum and trialNum
files=dir('*video.avi');
test=isempty(files);
if test==0; %you have already saved trials in this directory
    prefix=strcat(datestr(date, 'yymmdd'));
    expAccum=[];
    trialAccum=[];
    for i = 1:size(files,1),
        fileName = files(i,1).name;
        tok = regexp(fileName,[prefix '_(\d+)_(\d+)_video.avi'],'tokens');
        exp = tok{1}{1};
        trial = tok{1}{2};
        expAccum = [expAccum str2num(exp)];
        trialAccum = [trialAccum str2num(trial)];
    end;
    lastExp = max(expAccum);
    expNum=lastExp;
    lastTrial = max(trialAccum);
    trialNum=lastTrial+1;
else expNum=1; trialNum=1;
end

%initialize DAQ
numdatapts = sampleRate*numSecs;
s = daq.createSession('ni');

s.Rate=sampleRate;
s.addAnalogOutputChannel('Dev1','ao0','voltage'); %drives the red light
s.addAnalogInputChannel( 'Dev1','ai0','voltage'); %receives the same command which drives the red light
delay=zeros(((numdatapts/2)-1),1);
after=zeros(numafter,1);
t=(1:1:(numdatapts/2));
odoron(t)=5;
output1=vertcat(delay(1:end-1),odoron',after,0);
%output1=(output1)';
outputData(:,1)=output1;

delta = 1/3;
outputData1 = zeros(8000,1);
for i = 1:10
    outputData1 = [outputData1;delta*i*ones(8000,1)];
end
outputData1 = [outputData1;zeros(1000,1)];
s.queueOutputData(outputData1);
s.NotifyWhenDataAvailableExceeds =(length(outputData1)-100);
display('Data acquisition starting');
tic
%sampleRate*numSecs;

%start the listener in the background
dataListener = s.addlistener('DataAvailable', @saveData_optogenetics);
s.startBackground(); %start background is code to start data acquisition



