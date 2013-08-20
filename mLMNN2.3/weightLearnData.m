% This program is used to process training data for weight learning.
load('weightPerm.mat');
load('multiExpert.mat');

load('data/viper_features.mat');
load('finalFeature_5.mat');
%load('multiExpert.mat')
newFeature = [ux;feature];

training = {};
%a = [1 2 4 5 8 9 10 12 14 15 16 17 18 19 22 24 26 27 28 29];
a = [2 14 16 17 22 24 26 27 28];
selectSet = selectSet(a);
for i = 1:size(selectSet,2)
	myselect = selectSet{i}.rand;
	selectFeature = [];
	for j = 1:24
		st = myselect(j);
		selectFeature = [selectFeature;newFeature(st:(st+50),:)];
	end
	sample.training1 = selectFeature(:,perm(317:632));
	sample.training2 = selectFeature(:,perm(317:632)+632);
	training{i} = sample;
	clear sample;
end

save('weightLearn.mat','training','selectSet');