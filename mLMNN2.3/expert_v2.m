%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%        Author: Bi Sai                 %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%     Date: June 27th, 2013             %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Usage: Multi Distance Learning Metric %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load data
echo off;
clear all;
clc;
rand('seed',1);
setpaths

fprintf('Loading data ...\n');
%load('data/viper_features.mat'); 
load('expertData.mat');

% %random permuation to get sample training data
% permuation = randperm(632);
% sample1 = ux(:,permuation(1:316));
% sample2 = ux(:,permuation(1:316) + 632);
% 
% %save test data
% testData1 = ux(:,permuation(316:632));
% testData2 = ux(:,permuation(316:632) + 316);
% save('testData.mat', 'testData1', 'testData2');
% 
% sample = [sample1 sample2];
% id = 1:632;
% 
% %divide the positive pairs of training data into three parts
% training1 = [sample1(:,1:105) sample2(:,1:105)];
% label1 = [linspace(1,105,105) linspace(1,105,105)];
% id1 = [1:105 [1:105] + 316];
% 
% training2 = [sample1(:,106:210) sample2(:,1:210)];
% label2 = [linspace(106,210,105) linspace(106,210,105)];
% id2 = [106:210 [106:210]+316];
% 
% training3 = [sample1(:,211:end) sample2(:,211:end)];
% label3 = [linspace(211,316,106) linspace(211,316,106)];
% id3 = [211:316 [211:316] + 316];


%doing a LMNN to find the inital value for L
% [Ldim,Det]=lmnn2(sample,[1:316 1:316],1,'quiet',1,'maxiter',500,'checkup',0);
load('expertLdim.mat');

metric1 = Ldim;
metric2 = Ldim;
metric3 = Ldim;

%divide the negative pairs into three parts
% [x,y] = meshgrid(1:316, 317:632);
% negativePair = [x(:) y(:)];         %generate all pairs
% negativePair(1:316:(316*316)) = []; %remove positive pair
% 
% pairPerm = randperm(95440); 
% pairDist = zeros(632,632); %used to record the distance between each pair
% %random partition
% 
% for i = 1:33180
% 	point1 = sample(:,negativePair(pairPerm(i),1));
% 	point2 = sample(:,negativePair(pairPerm(i),2));
% 	pairDist(negativePair(pairPerm(i),1), negativePair(pairPerm(i),2)) = sum((metric1*(point1 - point2)).^2);
% end
% 
% for i = 33181:66360
% 	point1 = sample(:,negativePair(pairPerm(i),1));
% 	point2 = sample(:,negativePair(pairPerm(i),2));
% 	pairDist(negativePair(pairPerm(i),1), negativePair(pairPerm(i),2)) = sum((metric2*(point1 - point2)).^2);
% end
% 
% for i = 66361:95440
% 	point1 = sample(:,negativePair(pairPerm(i),1));
% 	point2 = sample(:,negativePair(pairPerm(i),2));
% 	pairDist(negativePair(pairPerm(i),1), negativePair(pairPerm(i),2)) = sum((metric3*(point1 - point2)).^2);
% end

load('negPair.mat');

%******************************************************************************
%Using LMNN
negGroup1 = [];
negGroup2 = [];
negGroup3 = [];

for iter = 1 : 5
    
 	metric1 = myLMNN(metric1,pairDist,id1,training1,label1,300,sample);
 	metric2 = myLMNN(metric2,pairDist,id2,training2,label2,300,sample);
 	metric3 = myLMNN(metric3,pairDist,id3,training3,label3,300,sample);
    %load('expertMetric.mat');
	%---------------------------------------------------------------------------------------
	%---------------------------------------------------------------------------------------
	%divide the positve pair into three parts again
	dist1 = sum((metric1 * (sample1 - sample2)).^2);
	dist2 = sum((metric2 * (sample1 - sample2)).^2);
	dist3 = sum((metric3 * (sample1 - sample2)).^2);

	dis = [dist1' dist2' dist3'];
	%find the minimum of each row
	%index records the index of minimum on each row
	[minimum, minIndex] = min(dis,[],2);  
    
    fprintf('Find group 1++++++++++++++++++++++++++\n');
	%find new group1
	index1 = find(minIndex == 1);
	training1 = [sample1(:,index1) sample2(:,index1)];
	id1 = [index1 index1 + 316];
	%find new group2
	
    fprintf('Find group 2++++++++++++++++++++++++++\n');
	index2 = find(minIndex == 2);
	training2 = [sample1(:,index2) sample2(:, index2)];
	id2 = [index2 index2 + 316];

    fprintf('Find group 3++++++++++++++++++++++++++\n');
	%find new group3
	index3 = find(minIndex == 3);
	training3 = [sample1(:,index3) sample2(:,index3)];
	id3 = [index3 index3 + 316];
    %---------------------------------------------------------------------------------------
	%---------------------------------------------------------------------------------------
	%divide negative pairs into three parts again
    
    
	negPair1 = sample(:,negativePair(:,1));
	negPair2 = sample(:,negativePair(:,2));
    fprintf('computing dist1+++++++++++++++++++++++++++++++++\n');
	dist1 = sum((metric1 *(negPair1 - negPair2)).^2);
    fprintf('computing dist2+++++++++++++++++++++++++++++++++\n');
	dist2 = sum((metric2 *(negPair1 - negPair2)).^2);
    fprintf('computing dist3+++++++++++++++++++++++++++++++++\n');
	dist3 = sum((metric3 *(negPair1 - negPair2)).^2);

	dis = [dist1' dist2' dist3'];
	%find maximum of each row
	[maximum, maxIndex] = max(dis,[],2);

    fprintf('Find negative paris of group 1++++++++++++++++++++++++++ \n');
	%find negative pairs that belong to group1
	index1 = find(maxIndex == 1);
	for i = 1 : size(index1,1)
		point1 = sample(:,negativePair(index1(i),1));
	    point2 = sample(:,negativePair(index1(i),2));
	    pairDist(negativePair(index1(i),1), negativePair(index1(i),2)) = sum((metric1*(point1 - point2)).^2);
	end
	negGroup1 = [sample(:,negativePair(index1,1)); sample(:,negativePair(index1,2))];
	
    
    fprintf('Find negative paris of group 2++++++++++++++++++++++++++  \n');
	%find negative pairs that belong to group2
	index2 = find(maxIndex == 2);
	for i = 1 : size(index2,1)
		point1 = sample(:,negativePair(index2(i),1));
	    point2 = sample(:,negativePair(index2(i),2));
	    pairDist(negativePair(index2(i),1), negativePair(index2(i),2)) = sum((metric2*(point1 - point2)).^2);
	end
    negGroup2 = [sample(:,negativePair(index2,1)); sample(:,negativePair(index2,2))];
	
    
    fprintf('Find negative paris of group 3+++++++++++++++++++++++++ \n');
	%find negative pairs that belong to group3
	index3 = find(maxIndex == 3);
	for i = 1 : size(index3,1)
		point1 = sample(:,negativePair(index3(i),1));
	    point2 = sample(:,negativePair(index3(i),2));
	    pairDist(negativePair(index3(i),1), negativePair(index3(i),2)) = sum((metric3*(point1 - point2)).^2);
	end
	negGroup3 = [sample(:,negativePair(index3,1)); sample(:,negativePair(index3,2))];
end

save('expertLMNN.mat', 'training1','training2','training3'...
    ,'negGroup1','negGroup2','negGroup3');


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Using SVM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:3
	%data points formed by positive pairs
	svmTraining1 = [training1(:,1:size(training1,2)/2);
	                training1(:,size(training1,2)/2 + 1 : end)];
	svmTraining2 = [training2(:,1:size(training2,2)/2);
	                training2(:,size(training3,2)/2 + 1 : end)];
	svmTraining3 = [training3(:,1:size(training3,2)/2);
	                training3(:,size(training3,2)/2 + 1 : end)];
    %data points formed by negative pairs
    svmTraining4 = [negGroup1 negGroup2 negGroup3]';

	svmTraining = [svmTraining1 svmTraining2 svmTraining3]';
	svmLabel = [ones(size(training1,2)/2,1);
	            ones(size(training2,2)/2,1) * 2;
	            ones(size(training3,2)/2,1) * 3;
	            ones(size(negGroup1,2),1);
	            ones(size(negGroup2,2),1) * 2;
	            ones(size(negGroup3,2),1) * 3];

	%parameter selection for positive pair 
    bestcv = 0;
    cv = {};
	for log2c = -1:3,
	  for log2g = -4:1,
	    cmd = ['-v 5 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
	    cv = svmtrain(svmLabel, svmTraining, cmd);
	    if (cv >= bestcv)
	      bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
	    end
	    fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
	  end
	end

    %do predication to the positive pairs and negative pairs
    predictData1 = [sample1;sample2];
    predictData2 = [sample(:,negativePair(:,1)); sample(:,negativePair(:,2))];
    predictData = [predictData1 predictData2]';

    [label, Accu, Prob] = svmpredict(svmLabel,predictData,cv, '-b',1);
  
    posLabel = label(1:316);
    %find group 1,only the first 316 points with label 1 belong to training1
    index1 = find(posLabel == 1);
    training1 = [sample1(:,index1) sample2(:,index1)];
    id1 = [index1  index1 + 316];

    %find group 2, only the first 316 points with label 2 belong to training2
    index2 = find(posLabel == 2);
    training2 = [sample1(:,index2) sample2(:,index2)];
    id2 = [index2 index2 + 316];

    %find group 3, only the first 316 points with label 3 belong to training3
    index3 = find(posLabel == 3);
    training3 = [sample1(:,index3) sample2(:,index3)];
    id3 = [index3 index3 + 316];

    %divide negative pairs into 3 groups
    negLabel = label(317:end);
    %find group1
    index1 = find(negLabel == 1);
    for i = 1 : size(index1,1)
		point1 = sample(:,negativePair(index1(i),1));
	    point2 = sample(:,negativePair(index1(i),2));
	    pairDist(negativePair(index1(i),1), negativePair(index1(i),2)) = sum((metric1*(point1 - point2)).^2);
	end

    %find group2
	index2 = find(negLabel == 2);
    for i = 1 : size(index2,1)
		point1 = sample(:,negativePair(index2(i),1));
	    point2 = sample(:,negativePair(index2(i),2));
	    pairDist(negativePair(index2(i),1), negativePair(index2(i),2)) = sum((metric2*(point1 - point2)).^2);
	end

	%find group3
	index3 = find(negLabel == 3);
    for i = 1 : size(index3,1)
		point1 = sample(:,negativePair(index3(i),1));
	    point2 = sample(:,negativePair(index3(i),2));
	    pairDist(negativePair(index3(i),1), negativePair(index3(i),2)) = sum((metric3*(point1 - point2)).^2);
	end



	%===============================================================================
	%===============================================================================
	%using LMNN again: code copied from above

	metric1 = myLMNN(metric1,pairDist,id1,training1,label1,500,sample); 
	metric2 = myLMNN(metric2,pairDist,id2,training2,label2,500,sample);
	metric3 = myLMNN(metric3,pairDist,id3,training3,label3,500,sample);

	%---------------------------------------------------------------------------------------
	%---------------------------------------------------------------------------------------
	%divide the positve pair into three parts again
	dist1 = sum((metric1 * (sample1 - sample2)).^2);
	dist2 = sum((metric2 * (sample1 - sample2)).^2);
	dist3 = sum((metric3 * (sample1 - sample2)).^2);

	dis = [dist1 dist2 dist3];
	%find the minimum of each row
	%index records the index of minimum on each row
	[minimum, minIndex] = min(dis,[],2);  

	%find new group1
	index1 = find(minIndex == 1);
	training1 = [sample1(:,index1) sample2(:,index1)];
	id1 = [index1 index1 + 316];
	%find new group2
	
	index2 = find(minIndex == 2);
	training2 = [sample1(:,index2) sample2(:, index2)];
	id2 = [index2 index2 + 316];

	%find new group3
	index3 = find(minIndex == 3);
	training3 = [sample1(:,index3) sample2(:,index3)];
	id3 = [index3 index3 + 316];
    %---------------------------------------------------------------------------------------
	%---------------------------------------------------------------------------------------
	%divide negative pairs into three parts again
	negPair1 = sample(:,negativePair(:,1));
	negPair2 = sample(:,negativePair(:,2));
	dist1 = sum((metric1 *(negPair1 - negPair2)).^2);
	dist2 = sum((metric2 *(negPair1 - negPair2)).^2);
	dist3 = sum((metric3 *(negPair1 - negPair2)).^2);

	dis = [dist1 dist2 dist3];
	%find maximum of each row
	[maximum, maxIndex] = max(dis,[],2);

	%find negative pairs that belong to group1
	index1 = find(maxIndex == 1);
	for i = 1 : size(index1,1)
		point1 = sample(:,negativePair(index1(i),1));
	    point2 = sample(:,negativePair(index1(i),2));
	    pairDist(negativePair(index1(i),1), negativePair(index1(i),2)) = sum((metric1*(point1 - point2)).^2);
	end
	negGroup1 = [sample(:,negativePair(index1,1)); sample(:,negativePair(index1,2))];
	
	%find negative pairs that belong to group2
	index2 = find(maxIndex == 2);
	for i = 1 : size(index2,1)
		point1 = sample(:,negativePair(index2(i),1));
	    point2 = sample(:,negativePair(index2(i),2));
	    pairDist(negativePair(index2(i),1), negativePair(index2(i),2)) = sum((metric2*(point1 - point2)).^2);
	end
    negGroup2 = [sample(:,negativePair(index2,1)); sample(:,negativePair(index2,2))];
	
	%find negative pairs that belong to group3
	index3 = find(maxIndex == 3);
	for i = 1 : size(index3,1)
		point1 = sample(:,negativePair(index3(i),1));
	    point2 = sample(:,negativePair(index3(i),2));
	    pairDist(negativePair(index3(i),1), negativePair(index3(i),2)) = sum((metric3*(point1 - point2)).^2);
	end
	negGroup3 = [sample(:,negativePair(index3,1)); sample(:,negativePair(index3,2))];

	%===============================================================================
	%===============================================================================
end
%}

%finished, got final metric
save('metric.mat', 'metric1', 'metric2', 'metric3');


















