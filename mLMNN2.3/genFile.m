% This program is used to write data from viper_feature.mat
% into files. The file is in SVM-LIGHT format. And for each
% file, the target value is a dimension in the feature Ldim
% * ux. And the feature pairs are the original features.

load('uxLdim2.mat');
load('data/viper_features.mat');
load('perm.mat');
perm = permuation;

feature = ux;
metric = uxLdim;
% training set
training1 = feature(:,perm(1:270));
training2 = feature(:,perm(1:270)+632);
% validation set
validation1 = feature(:,perm(271:316));
validation2 = feature(:,perm(271:316)+632);
% test set
test1 = feature(:,perm(317:632));
test2 = feature(:,perm(317:632)+632);

mTraining1 = metric * training1;
mTraining2 = metric * training2;
mValidation1 = metric * validation1;
mValidation2 = metric * validation2;
mTest1 = metric * test1;
mTest2 = metric * test2; 

for dim = 1:size(mTraining1,1)
	str = int2str(dim);

	fileName = strcat('GBRT_DATA/dimTrain',str,'.gbrt');
	fileID = fopen(fileName,'w+');
	writeFile(fileID,[training1 training2],dim,[mTraining1 mTraining2]);
	fclose(fileID);

	fileName = strcat('GBRT_DATA/dimValidation',str,'.gbrt');
	fileID = fopen(fileName,'w+');
	writeFile(fileID,[validation1 validation2], dim, [mValidation1 mValidation2]);
	fclose(fileID);

	fileName = strcat('GBRT_DATA/dimTest',str,'.gbrt');
	fileID = fopen(fileName,'w+');
	writeFile(fileID,[test1 test2], dim, [mTest1 mTest2]);
	fclose(fileID);	
end


save('GBRT_DATA/data.mat','training1','training2','test1','test2','validation1','validation2');



















