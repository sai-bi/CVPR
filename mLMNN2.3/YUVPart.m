%This program is used to do LMNN with YUV part
load('finalFeature.mat');
load('permuation.mat');

sample1 = YUVpart(:,permuation(1:316));
sample2 = YUVpart(:,permuation(1:316) + 632);
YUVSample = [sample1 sample2];
YUVTest1 = YUVpart(:,permuation(317:632));
YUVTest2 = YUVpart(:,permuation(317:632) + 632);

[LdimYUV,det1] = multiLMNN(YUVSample,[1:316 1:316],1,'maxiter',600,'checkup',0,'outdim',80);
save('YUVPart.mat','YUVTest1','LdimYUV','YUVTest2','YUVSample');
