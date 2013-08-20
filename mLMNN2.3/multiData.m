load('finalFeature_5.mat');
load('data/viper_features.mat');

newFeature = [ux; feature];


%save('randSelect.mat','training1','training2','test1','test2');

selectSet = {};
perm = randperm(632);
for i = 1:30
    select = randsample(1221,24);

    %selectFeature = zeros(600,1264);
    selectFeature = [];
    for j = 1:24
        st = select(j);
        selectFeature = [selectFeature;newFeature(st:(st+50),:)];
    end

    sample.rand = select;
    sample.training1 = selectFeature(:,perm(1:316));
    sample.training2 = selectFeature(:,perm(1:316)+632);
    sample.test1 = selectFeature(:,perm(317:632));
    sample.test2 = selectFeature(:,perm(317:632)+632);
    selectSet{i} = sample;
    clear sample;
end

save('randSelect.mat','selectSet','perm');
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
