function [newState] = BorderChoice(self,state,f,df,baseline,fs,originalFs)
m = self.model.BorderChoice;
newStateProb = m.during(df,f).*originalFs./fs;
newStateProb(df ==0 & abs(f-baseline)<0.001) = 0;
newStateProb(~state) = 0;

newState = rand(size(f))<newStateProb;

end