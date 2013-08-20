echo off;
clear all;
clc;
rand('seed',1);
setpaths

load('testData.mat');
% fprintf('Loading data ...\n');
% load('data/viper_features.mat');
% load('test.mat');
% % permuation = randperm(632);
% % sample1 = ux(:,permuation(1:316));
% % sample2 = ux(:,permuation(1:316) + 632);
% % 
% % %save test data
% % testData1 = ux(:,permuation(317:632));
% % testData2 = ux(:,permuation(317:632) + 632);
% % save('testData.mat', 'testData1', 'testData2');
% % 
% % sample = [sample1 sample2];
% % id = 1:632;

[Ldim2,Det]= multiLMNN(sample,[1:316 1:316],1,'maxiter',450,'checkup',0,'outdim',60);

save('Ldim2.mat','Ldim2');






