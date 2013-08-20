%random permuation to get sample training data
load('data/viper_features.mat');
permuation = randperm(632);
sample1 = ux(:,permuation(1:316));
sample2 = ux(:,permuation(1:316) + 632);

%save test data
testData1 = ux(:,permuation(317:632));
testData2 = ux(:,permuation(317:632) + 632);


sample = [sample1 sample2];
id = 1:632;

%divide the positive pairs of training data into three parts
training1 = [sample1(:,1:105) sample2(:,1:105)];
label1 = [linspace(1,105,105) linspace(1,105,105)];
id1 = [1:105 [1:105] + 316];

training2 = [sample1(:,106:210) sample2(:,106:210)];
label2 = [linspace(106,210,105) linspace(106,210,105)];
id2 = [106:210 [106:210]+316];

training3 = [sample1(:,211:end) sample2(:,211:end)];
label3 = [linspace(211,316,106) linspace(211,316,106)];
id3 = [211:316 [211:316] + 316];

save('expertData.mat','sample1','sample2','sample','testData1','testData2'...
     ,'id','id1','id2','id3','training1','training2','training3','label1','label2'...
     ,'label3');
 
save('expertSample.mat','sample');