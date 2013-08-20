% This program is used to learn weight for metrics.

load('weightLearn.mat');

trainNum = size(training{1}.training1,2);
metricNum = size(selectSet,2);
sample1 = [];
sample2 = [];
% for i = 1: trainNum
 	trainIns1 = [];
 	trainIns2 = [];
	for j = 1:metricNum
        fprintf('j: %d, i: %d \n',j, i);
		temp1 = (selectSet{j}.metric) * training{j}.training1;
		trainIns1 = [trainIns1; temp1];
		temp2 = (selectSet{j}.metric) * training{j}.training2;
		trainIns2 = [trainIns2; temp2];
	end
	sample1 = [sample1 trainIns1];
	sample2 = [sample2 trainIns2];
%end

[outputMetric, Det]= lmnn2([sample1 sample2],[1:trainNum 1:trainNum],1,'maxiter',15000,'checkup',0,'outdim',80);


save('ensemble.mat','outputMetric');





















