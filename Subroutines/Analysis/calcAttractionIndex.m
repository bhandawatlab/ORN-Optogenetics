function attractNdx = calcAttractionIndex(self,lab,border)
% calcAttractionIndex  calculates the attraction index
%
%   Inputs: self = fly object
%           lab = label ('' = body, 'H' = head)
%           border = light arena border (in cm)
%
%   Output: attractNdx = structure for attraction index and raw amount of
%       time spent inside and outside
%   

fe = getFirstEntry(self,lab,border);
for i = 1:self.nFly
    totInside = sum(self.(['r' lab])(i,fe(i):end)<border);
    totOutside = sum(self.(['r' lab])(i,fe(i):end)>=border);
    attractNdx.during(i) = totInside./(totInside+totOutside);
    attractNdx.duringIn(i) = totInside;
    attractNdx.duringTot(i) = totInside+totOutside;
    
    totInside = sum(self.(['r' lab])(i,1:fe(i))<border);
    totOutside = sum(self.(['r' lab])(i,1:fe(i))>=border);
    attractNdx.before(i) = totInside./(totInside+totOutside);
    attractNdx.beforeIn(i) = totInside;
    attractNdx.beforeTot(i) = totInside+totOutside;
end

end