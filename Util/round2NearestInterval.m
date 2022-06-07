function [out] = round2NearestInterval(A,interval)
if isempty(interval)
    interval = 5;
end
numPow = floor(log10(abs(A)));
out = ceil(A./((10.^(numPow-1)).*interval)).*((10.^(numPow-1)).*interval);
end