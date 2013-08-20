load('finalFeature.mat');
load('permuation.mat');

sample1 = HSVpart(:,permuation(1:316));
sample2 = HSVpart(:,permuation(1:316)+632);
HSVsample = [sample1 sample2];
HSVTest1 = HSVpart(:,permuation(317:632));
HSVTest2 = HSVpart(:,permuation(317:632)+632);

[LdimHSV,det1] = multiLMNN(HSVsample,[1:316 1:316],1,'maxiter',600,'checkup',0,'outdim',80);

save('HSVPart.mat','LdimHSV','HSVTest1','HSVTest2','HSVsample');
