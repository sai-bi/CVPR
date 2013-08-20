%This program is used to test myLMNN.myLMNN

echo off;
clear all;
clc;
rand('seed',1);
setpaths

fprintf('Loading data ...\n');
load('test.mat');
load('pair.mat');

% pairDist = zeros(316,316);
% for i = 1 : 316
% 	for j = 317 : 632
% 		temp1 = Ldim * (sample(:,i) - sample(:,j));
% 		pairDist(i,j - 316) = temp1' * temp1;
% 	end
% end

% save('pair.mat','pairDist');

% load('test.mat');

metric = myLMNN(eye(631), pairDist, 1:632,sample,[1:316 1:316],300,sample);





