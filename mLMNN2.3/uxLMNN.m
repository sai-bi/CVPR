%load('permuation.mat');
load('data/viper_features.mat');

permuation = randperm(632);
sample1 = ux(:,permuation(1:316));
sample2 = ux(:,permuation(1:316)+632);

uxTest1 = ux(:,permuation(317:632));
uxTest2 = ux(:,permuation(317:632)+632);

save('uxData.mat','sample1','sample2','uxTest1','uxTest2');
save('perm.mat','permuation');