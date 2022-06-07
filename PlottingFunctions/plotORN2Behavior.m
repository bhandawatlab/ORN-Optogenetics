function [mat] = plotORN2Behavior(labCov2Cons,allData)

nGen2cons = numel(labCov2Cons);
nBehaviors = numel(cell2mat(allData{1, 4}.values(:,1)'));
behaviorLabel = cell(3,nBehaviors);
mat = zeros(nBehaviors,numel(labCov2Cons));
for i = 1:numel(allData)
    for j = 1:numel(allData{i}.labCov)
        ORNNdx = strcmpi(allData{i}.labCov{j},labCov2Cons);
        if sum(ORNNdx)>0
            currKin = cell2mat(allData{i}.values(:,j)');
            mat(:,ORNNdx) = currKin;
        end
    end
    k = 1;
    for stateType = 1:size(allData{i}.stateType,2)
        for ROI = 1:numel(allData{i}.ROI)
            behaviorLabel(:,k) = [allData{i}.stateType(:,stateType);...
                allData{i}.ROI(:,ROI)];
            k = k+1;
        end
    end
end

behaviorLabel = strrep(behaviorLabel,'sharp turns','ST');
behaviorLabel = strrep(behaviorLabel,'curved walks','CW');

thresh = 0.15;
genPos = linspace(nGen2cons,nBehaviors-nGen2cons,nGen2cons);
figure;set(gcf,'Position',[2 42 838 924])
scatter(2.*ones(1,nGen2cons),genPos,'ok');hold on;
scatter(12.*ones(1,nBehaviors),1:nBehaviors,'ok');
text(0.5.*ones(1,nGen2cons),genPos,labCov2Cons)
text(12.5.*ones(1,nBehaviors),1:nBehaviors,join(behaviorLabel'))
for gen = 1:numel(labCov2Cons)
    for beh = 1:nBehaviors
        if mat(beh,gen)>(1+thresh)
            plot([2,12],[genPos(gen),beh],'r')
        elseif mat(beh,gen)<(1-thresh)
            plot([2,12],[genPos(gen),beh],'b')
        end
    end
end
xlim([0 15])

mat2 = mat;
mat2(mat>(1+thresh)) = 1;
mat2(mat<(1-thresh)) = -1;
mat2(mat>=(1-thresh) & mat<=(1+thresh)) = 0;

figure;set(gcf,'Position',[2 42 838 924])
imagesc(mat2)
colormap([0.8500 0.3250 0.0980;1,1,1;0.3010 0.7450 0.9330]);
t2 = num2str(round(mat(:),2));
[X,Y] = meshgrid(1:nGen2cons,1:nBehaviors);
text(X(:),Y(:),t2, 'HorizontalAlignment', 'Center','Color',[0 0 0])
set(gca,'ytick',1:nBehaviors,'yticklabel',join(behaviorLabel'));
set(gca,'xtick',1:nGen2cons,'xticklabel',labCov2Cons);

end