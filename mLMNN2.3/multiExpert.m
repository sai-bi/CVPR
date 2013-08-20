load('randSelect.mat');

for i = 1:size(selectSet,2)
	fprintf('%d+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n',i);
	fprintf('%d+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n',i);
	fprintf('%d+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n',i);

	training = [selectSet{i}.training1 selectSet{i}.training2];
	label = [1:316 1:316];
	[metric1,Det1]= lmnn2(training,label,1,'maxiter',15000,'checkup',0,'outdim',100);
	selectSet{i}.metric = metric1;
end

save('multiExpert.mat','selectSet');