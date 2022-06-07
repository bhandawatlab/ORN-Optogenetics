function [r] = sampleKinematics_Boundary(self,state,n)
m = self.model.boundary;
nKin = 2;

if state == 1% before
    muAll = repmat(m.baseline(:,1)',n,1);
    sigmaAll = repmat(m.baseline(:,2)',n,1);
    qBounds = repmat(m.qBoundsBaseline(2,:),n,1);
elseif state == 2
    muAll = repmat(m.during(:,1)',n,1);
    sigmaAll = repmat(m.during(:,2)',n,1);
    qBounds = repmat(m.qBounds(2,:),n,1);
end

r = lognrnd(muAll,sigmaAll);
r(r>qBounds) = qBounds(r>qBounds);

end