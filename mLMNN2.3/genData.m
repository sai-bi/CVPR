echo off;
clear all;
clc;
rand('seed',1);
setpaths

fprintf('Loading data ...\n');
load('finalFeature.mat');

%feature = mixFeature;
permuation = randperm(632);
save('permuation.mat','permuation');

sample1 = feature(:,permuation(1:316));
sample2 = feature(:,permuation(1:316) + 632);

%save test data
testData1 = feature(:,permuation(317:632));
testData2 = feature(:,permuation(317:632) + 632);
% save('testData.mat', 'testData1', 'testData2');

sample = [sample1 sample2];
id = 1:632;

save('testData.mat','sample','testData1','testData2');