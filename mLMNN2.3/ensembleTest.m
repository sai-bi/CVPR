% This program is used to test the ensemble metric.

load('weightTest2.mat');
load('ensemble.mat');

metricNum = size(selectSet,2);
testNum = size(test{1}.test1,2);
test1 = [];
test2 = [];
for i = 1:metricNum
	temp1 = selectSet{i}.metric * test{i}.test1;
	temp2 = selectSet{i}.metric * test{i}.test2;

	test1 = [test1;temp1];
	test2 = [test2;temp2];
end

%test1 = outputMetric * test1;
%test2 = outputMetric * test2;

x = 1:testNum;
y2 = zeros(1,testNum);
for i = 1 : testNum
	temp2 = bsxfun(@minus, test2,test1(:,i));
    temp2 = sum(temp2.^2,1);
	temp1 = bsxfun(@minus, test1,test2(:,i));
    temp1 = sum(temp1.^2,1);
	temp = temp1 + temp2;
	[mininum, minIndex] = sort(temp);
    minIndex = find(minIndex == i);
    %display(minIndex);
	y2 = y2 + [zeros(1,minIndex-1) ones(1,testNum + 1-minIndex)];
end

y2 = y2/testNum * 100;

Z = trapz(x,y2);
fprintf('%f\n',Z/((testNum - 1)*100));
plot(x,y2,'color','g');
hold on;

















