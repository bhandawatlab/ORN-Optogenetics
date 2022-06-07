function [r] = sampleKinematicsJoint_Boundary(self,state,n)
m = self.model.boundary;
nKin = 2;

if state == 1% before
    muAll = repmat(m.baseline(:,1)',n,1);
    sigmaAll = repmat(m.baseline(:,2)',n,1);
    
    corrMat = repmat(m.corrBaseline,1,1,n);
    qBound = repmat(m.qBoundsBaseline,1,1,n);
elseif state == 2
    muAll = repmat(m.during(:,1)',n,1);
    sigmaAll = repmat(m.during(:,2)',n,1);
    
    corrMat = repmat(m.corr,1,1,n);
    qBound = repmat(m.qBounds,1,1,n);
end

% with sigma
%--------------------------------------------------------------------------
r = zeros(n,nKin);
template = repmat([1 -1],nKin,1);
for i = 1:n
    try 
        y = MvLogNRand( muAll(i,:) , sigmaAll(i,:)' , 1 , corrMat(:,:,i) );
    catch
        y = MvLogNRand( muAll(i,:) , sigmaAll(i,:)' , 1 , abs(corrMat(:,:,i)) );
    end
    inBound = all(sign((y'-qBound(:,:,i)))==template,'all');
    while ~inBound
        try 
            y = MvLogNRand( muAll(i,:) , sigmaAll(i,:)' , 1 , corrMat(:,:,i) );
        catch
            y = MvLogNRand( muAll(i,:) , sigmaAll(i,:)' , 1 , abs(corrMat(:,:,i)) );
        end
        inBound = all(sign((y'-qBound(:,:,i)))==template,'all');
    end
    r(i,:) = y;
end
%r = exp(mvnrnd(muAll,sigma));

% % without sigma
% %--------------------------------------------------------------------------
% r = exp(muAll);

end