function [out] = getFilteredResponse(I,nFly,b)
n = numel(b)-1;

I(I<0) = 0;
I_pad = [zeros(nFly,n) I];%I(:,3:end);
out = conv2(I_pad,flip(b(2:end),1)','valid')+b(1);
out(out<0) = 0;

end
