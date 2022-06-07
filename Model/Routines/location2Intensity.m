function [V_out,I_out] = location2Intensity(rPos,x,I,convIV)

[~,idxB] = min(abs(rPos-x),[],2);
I_out =  I(idxB);
V_out = convIV(I_out);
V_out(V_out==convIV(0)) = 0;

end