function [spk] = getSpikeRate(I,nFly,bLFP,bRate)
n1 = 500;n2 = 500;
%l = size(I,2);
totN = n1+n2;
%rat = fs./filterSamp;

I(I<0) = 0;
I_pad = [zeros(nFly,totN-2) I];%I(:,3:end);
LFP = conv2(I_pad,flip(bLFP(2:end),1)','valid')+bLFP(1);%convolve2
spk = conv2(LFP,flip(bRate(2:end),1)','valid')+bRate(1);%convolve2
spk(spk<0) = 0;


% %I = [zeros(nFly,ceil(totN.*rat)-l), I];
% %I = [zeros(nFly,ceil(totN.*rat)-1), I];
% if l>1
%     I = resample(I',filterSamp,fs)';
% end
% %I = smoothdata(I,2,'lowess',15);% maybe use 10?
% I(I<0) = 0;
% I_pad = [zeros(nFly,totN-2) I];%I(:,3:end);
% 
% LFP = conv2(I_pad,flip(bLFP(2:end),1)','valid')+bLFP(1);%convolve2
% spk = conv2(LFP,flip(bRate(2:end),1)','valid')+bRate(1);%convolve2
% spk(spk<0) = 0;
% spk = resample(spk',fs,filterSamp)';

end

% function [spk] = getSpikeRate(I,nFly,bLFP,bRate,fs,filterSamp)
% n1 = 500;n2 = 500;
% l = size(I,2);
% totN = n1+n2;
% rat = fs./filterSamp;
% 
% I = [zeros(nFly,ceil(totN.*rat)-l), I];
% I = resample(I',filterSamp,fs)';
% %I = smoothdata(I,2,'lowess',15);% maybe use 10?
% I(I<0) = 0;
% 
% LFP = conv2(flip(bLFP(2:end),1)',I)+bLFP(1);%convolve2
% LFP = LFP(:,1:n1+n2);
% spk = conv2(flip(bRate(2:end),1)',LFP)+bRate(1);%convolve2
% spk = spk(:,n1:n1+n2);
% spk(spk<0) = 0;
% 
% end