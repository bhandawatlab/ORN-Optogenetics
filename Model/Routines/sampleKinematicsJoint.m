function [r] = sampleKinematicsJoint(self,state,t,f,df,tt)
%t(t>110) = 110;
t(t>10) = 10;
%t = t.*0+1;
n = numel(tt);

nKin = numel(state);

corrMat = zeros(nKin,nKin,n);
muAll = zeros(n,nKin);
sigmaAll = zeros(n,nKin);
qBound = zeros(nKin,2,n);
for s = 1:nKin
    m = self.model.params{state(s)};
    
    y = m.h(m.tt,m.params.y,t);
    phi = m.h(m.tt,m.params.phi,t);
    a1 = m.h(m.tt,m.params.a1,t);
    a2 = m.h(m.tt,m.params.a2,t);
    b1 = m.h(m.tt,m.params.b1,t);
    b2 = m.h(m.tt,m.params.b2,t);

    % mapping f/df to kinematic distribution
    %----------------------------------------------------------------------
    %mu = m.g(m.params.c1,m.params.c2,y,phi,a1,a2,b1,b2,f,df);
    mu = m.g2{1}(df,f);
    %----------------------------------------------------------------------
    mu(t<0) = m.baseline(1);% when t<0, it means before first entry
    mu(tt>0 & tt<=2) = m.baselineEarly(1);
    mu(tt>2) = m.baselineLate(1);
    
    %----------------------------------------------------------------------
    %Sigma = m.sigma(abs(mu),m.b);
    Sigma = m.g2_std{1}(df,f);
    %----------------------------------------------------------------------
    Sigma(t<0) = m.baseline(2);% when t<0, it means before first entry
    Sigma(tt>0 & tt<=2) = m.baselineEarly(2);
    Sigma(tt>2) = m.baselineLate(2);
    
    muAll(:,s) = mu;
    sigmaAll(:,s) = Sigma;%./10
    
    corrMat(:,:,t<0) = repmat(m.corrBaseline,1,1,sum(t<0));%before first entry
    corrMat(:,:,t>=0) = repmat(m.corr,1,1,sum(t>=0));%after first entry
    
    qBound(:,:,t<0) = repmat(m.qBoundsBaseline,1,1,sum(t<0));%before first entry
    qBound(:,:,t>=0) = repmat(m.qBounds,1,1,sum(t>=0));%before first entry
end

% with sigma
%--------------------------------------------------------------------------
r = zeros(n,nKin);
template = repmat([1 -1],nKin,1);
x = 1;
for i = 1:n
    try 
        y = MvLogNRand( muAll(i,:) , sigmaAll(i,:)' , 1 , corrMat(:,:,i) );
    catch
        try
            y = MvLogNRand( muAll(i,:) , sigmaAll(i,:)' , 1 , abs(corrMat(:,:,i)) );
        catch
            y = MvLogNRand( muAll(i,:) , sigmaAll(i,:)'./2 , 1 , abs(corrMat(:,:,i)) );
        end
    end
    inBound = all(sign((y'-qBound(:,:,i)))==template,'all');
    while ~inBound
        x = x+1;
        try 
            y = MvLogNRand( muAll(i,:) , sigmaAll(i,:)' , 1 , corrMat(:,:,i) );
        catch
            try
                y = MvLogNRand( muAll(i,:) , sigmaAll(i,:)' , 1 , abs(corrMat(:,:,i)) );
            catch
                y = MvLogNRand( muAll(i,:) , sigmaAll(i,:)'./2 , 1 , abs(corrMat(:,:,i)) );
            end
        end
        inBound = all(sign((y'-qBound(:,:,i)))==template,'all');
        if x>100
            y(y>qBound(:,2,i)') = qBound(y>qBound(:,2,i)',2,i);
            y(y<qBound(:,1,i)') = qBound(y<qBound(:,1,i)',2,i);
        end
    end
    r(i,:) = y;
end

% % without sigma
% %--------------------------------------------------------------------------
% r = exp(muAll);

end