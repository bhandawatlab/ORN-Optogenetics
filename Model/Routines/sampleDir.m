function [direction] = sampleDir(self,state,t,f,df)
%t = t.*0+1;
n = numel(t);

nKin = numel(state);

for s = 1:nKin
    m = self.model.params{state(s)};
    
    corrDir = m.dirCorr.corrFun(df,f);
    corrDir(t<0) = m.dirCorr.baselineCorr;
    
    direction = (rand(n,1)<=corrDir).*2-1; % set to -1 for diff dir and 1 for same dir
    
end


end