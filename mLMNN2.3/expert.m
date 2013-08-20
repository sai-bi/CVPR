%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%        Author: Bi Sai                 %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%     Date: June 23th, 2013             %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Topic: Multi Distance Learning Metric %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load data
echo off;
clear all;
clc;
rand('seed',1);
setpaths

fprintf('Loading data ...\n');
load('data/viper_features.mat');    

%First Part, generating data.
%The data is divided into three parts.
id = 1:632; %the id of each data point

perm = randperm(632);
sample1 = ux(perm(1:316));
sample2 = ux(perm(1:316) + 632);
sample1ID = ux(perm(1:316));
sample2ID = ux(perm(1:316) + 632);

%based on the id of the pair, we can identify their belonging
pairBelong = zeros(1264,1264);

%first group
training1part1 = ux(perm(:,1:105));
training1part2 = ux(perm(:,1:105) + 632);
training1 = [training1part1 training1part2];
label1 = [linspace(1:105:105) linspace(1:105:105)];
id1 = [id(perm(:,1:105)) id(perm(:,1:105) + 632)];

for i = 1:105
	pairBelong(id1(i,1), id1(i,2) + 632) = 1;
end

%second group
training2part1 = ux(perm(:,106:210));
training2part2 = ux(perm(:,106:210) + 632);
training2 = [training2part1 training2part2];
label2 = [linspace(106:210:105) linspace(106:210:105)];
id2 = [id(perm(:,106:210) id(perm(:,106:210) + 632)];

for i = 1:105
	pairBelong(id2(i,1), id2(i,2) + 632) = 2;
end

%third group
training3part1 = ux(perm(:,211:end));
training3part2 = ux(perm(:,211:end) + 632);
training3 = [training3part1 training3part2];
label3 = [linspace(211:316:106) linspace(211:316:106)];
id3 = [id(perm(:,211:end)) id(perm(:,211:end) + 632)];

for i = 1:106
	pairBelong(id3(i,1), id3(i,2)+ 632 ) = 3;
end


%Test
testpart1 = ux(perm(:,317:632));
testpart2 = ux(perm(:,317:632) + 632);
test = [testpart1 testpart2];
testlabel = [linspace(317,632,316) linspace(317,632,316)];


%partition the negative pairs
%there are C(316,2) - 316negative pairs in total
[x,y] = meshgrid(sample1ID, sample2ID);
negativePair = [x(:) y(:)];       %generate all pairs
negativePair(1:316:(316*316)) = []; %remove positive pair

pairPerm = randperm(1:95440);


for i = 1:33180
	pairBelong(negativePair(pairPerm(i),1), negativePair(pairPerm(i),2)) = 1;
end

for i = 33181:66360
	pairBelong(negativePair(pairPerm(i),1), negativePair(pairPerm(i),2)) = 2;
end

for i = 66361:95440
	pairBelong(negativePair(pairPerm(i),1), negativePair(pairPerm(i),2)) = 3;
end

negativeGroup1 = [];
negativeGroup2 = [];
negativeGroup3 = [];

%Main Loop 
for i = 1:10
	[Lmetric1,Det1]=lmnn2(pairBelong,training1,label1,1,'quiet',1,'maxiter',500,'checkup',0);
	[Lmetric2,Det2]=lmnn2(pairBelong,training2,label2,1,'quiet',1,'maxiter',500,'checkup',0);
	[Lmetric3,Det3]=lmnn2(pairBelong,training3,label3,1,'quiet',1,'maxiter',500,'checkup',0);

	%calcuate the distance between positive paris
	posPairDist1 = sum((Lmetric1 * (sample1 - sample2)).^2); 
	posPairDist2 = sum((Lmetric2 * (sample1 - sample2)).^2); 
	posPairDist3 = sum((Lmetric3 * (sample1 - sample2)).^2); 

	%find the smallest among these three distances
	posPairDist = [posPairDist1 posPairDist2 posPairDist3];
	[smallest, index] = min(posPairDist,[],2);
	
    posPairBelong = [sample1ID sample2ID];
	
	%positive pair that belongs to group 1
	index1 = find(index == 1);
	training1part1 = sample1(index1);
	training1part2 = sample2(index1);
	training1 = [training1part1 training1part2];
	id1 = [sample1ID(index1) sample2ID(index1)];
    
    pairBelong = [];
    for i = 1 : size(id1,1)
    	pairBelong(id1(i,1), id1(i,2)) = 1;
    end

    %positive pair that belongs to group 2
	index1 = find(index == 2);
	training2part1 = sample1(index1);
	training2part2 = sample2(index1);
	training2 = [training2part1 training2part2];
	id2 = [sample1ID(index1) sample2ID(index1)];
   
    for i = 1:size(id2,1)
    	pairBelong(id2(i,1), id2(i,2)) = 2;
    end

	%positive pair that belongs to group 3
	index1 = find(index == 3);
	training3part1 = sample1(index1);
	training3part2 = sample2(index1);
	training3 = [training3part1 training3part2];
	id3 = [sample1ID(index1)  sample2ID(index1)];

    %for pair 3 
	for i = 1:size(id3,1)
		pairBelong(id3(i,1), id3(i,2)) = 3;
	end

	%calculate the distance between negative pairs
	firstColumn = ux(:,negativePair(:,1));
	secondColumn = ux(:,negativePair(:,2));

	negPairDist1 = sum((Lmetric1 * (firstColumn - secondColumn)).^2);
	negPairDist2 = sum((Lmetric2 * (firstColumn - secondColumn)).^2);
	negPairDist3 = sum((Lmetric3 * (firstColumn - secondColumn)).^2);

	negPairDist = [negPairDist1 negPairDist2 negPairDist3];
	
	[largest, index] = max(negPairDist, [], 2);

	index1 = find(index == 1);
    for i = 1:size(index1,1)
		pairBelong(negativePair(index1(i),1), negativePair(index1(i),2)) = 1;
	end
	negativeGroup1 = [ux(:,negativePair(index1,1)); ux(:, negativePair(index1,2))];

	index1 = find(index == 2);
	for i = 1:size(index1,1)
		pairBelong(negativePair(index(i),1), negativePair(index2(i),2)) = 2;
	end
	negativeGroup2 = [ux(:,negativePair(index1,1)); ux(:, negativePair(index1,2))];

	index1 = find(index == 3);
	for i = 1:size(index1,1)
		pairBelong(negativePair(index(i),1), negativePair(index2(i),2)) = 3;
	end
	negativeGroup3 = [ux(:,negativePair(index1,1)); ux(:, negativePair(index1,2))];

end 


%Using SVM to do classification 
for i = 1:3 
	positive1 = [training1part1;training1part2];
	positive2 = [training2part1;training2part2];
	positive3 = [training3part1;training3part2];
	positive = [positive1 positive2 positive3];

	negative = [negativeGroup1 negativeGroup2 negativeGroup3];

	posLabel = [ones(1,size(positive1,2) ones(1,size(positive2,2)) * 2  ones(1,size(positive3,2)) * 2]';
    negLabel = [ones(1,size(negativeGroup1,2)) ones(1,size(negativeGroup2,2))*2 ones(1,size(negativeGroup3,3))*3]';

    %parameter selection for positive pair 
    bestcv = 0;
    cv = {};
	for log2c = -1:3,
	  for log2g = -4:1,
	    cmd = ['-v 5 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
	    cv = svmtrain(posLabel, positive', cmd);
	    if (cv >= bestcv)
	      bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
	    end
	    fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
	  end
	end

    [posLabel, posAccu, posProb] = svmpredict(posLabel, positive, cv, '-b',1);
    
    group1index = find(posLabel == 1);
    group1 = positive(:,group1index);
    training1part1 = group1(1:size(group1,1)/2, :);
    training1part2 = group1((size(group1,1)/2 + 1) : end, :);
    training1 = [training1part1 training1part2];
    label = linspace(1,size(training1part1,2), size(training1part1,2));
    label1 = [label label];

    group2index = find(posLabel == 2);
    group2 = positive(:,group2index);
    training2part1 = group2(1:size(group2,1)/2, :);
    training2part2 = group2((size(group2,1)/2 + 1) : end, :);
    training2 = [training2part1 training2part2];
    label = linspace(1,size(training1part1,2), size(training1part1,2)) + size(training1part1,2);
    label2 = [label label];

    group3index = find(posLabel == 3);
    group3 = positive(:,group3index);
    training3part1 = group3(1:size(group3,1)/2, :);
    training3part2 = group3((size(group3,1)/2 + 1) : end, :);
    training3 = [training3part1 training3part2]; 
    label = linspace(1,size(training3part1,2), size(training3part1,2)) + size(training1part1) + size(training2part1);
    label3 = [label label];

   

    %parameter selection for negative pair
    bestcv = 0;
    cv = {};
	for log2c = -1:3,
	  for log2g = -4:1,
	    cmd = ['-v 5 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
	    cv = svmtrain(negLabel, negative', cmd);
	    if (cv >= bestcv)
	      bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
	    end
	    fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
	  end
	end

	negative = [ux(negativePair(:,1)); ux(negativePair(:,2))];
	[negLabel, negAccu, negProb] = svmpredict(negLabel,negative,cv,'-b',1);

	group1index = find(negLabel == 1);
	for i = 1 : size(group1index,1)
		pairBelong(negativePair(i,1), negativePair(i,2)) = 1;
	end

    group2index = find(negLabel == 2);
    for i = 1 :  size(group2index,1)
    	pairBelong(negativePair(i,1), negativePair(i,2)) = 2;
    end

    group3index = find(label == 3);
    for i = 1 : size(group3index,1)
    	pairBelong(negativePair(i,1), negativePair(i,2)) = 3;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Using LMNN again:

    [Lmetric1,Det1]=lmnn2(pairBelong,training1,label1,1,'quiet',1,'maxiter',500,'checkup',0);
	[Lmetric2,Det2]=lmnn2(pairBelong,training2,label2,1,'quiet',1,'maxiter',500,'checkup',0);
	[Lmetric3,Det3]=lmnn2(pairBelong,training3,label3,1,'quiet',1,'maxiter',500,'checkup',0);

	%calcuate the distance between positive paris
	posPairDist1 = sum((Lmetric1 * (sample1 - sample2)).^2); 
	posPairDist2 = sum((Lmetric2 * (sample1 - sample2)).^2); 
	posPairDist3 = sum((Lmetric3 * (sample1 - sample2)).^2); 

	%find the smallest among these three distances
	posPairDist = [posPairDist1 posPairDist2 posPairDist3];
	[smallest, index] = min(posPairDist,[],2);
	
    posPairBelong = [sample1ID sample2ID];
	
	%positive pair that belongs to group 1
	index1 = find(index == 1);
	training1part1 = sample1(index1);
	training1part2 = sample2(index1);
	training1 = [training1part1 training1part2];
	id1 = [sample1ID(index1) sample2ID(index2)];
    
    pairBelong = [];
    for i = 1 : size(id1,1)
    	pairBelong(id1(i,1), id1(i,2)) = 1;
    end

    %positive pair that belongs to group 2
	index1 = find(index == 2);
	training2part1 = sample1(index1);
	training2part2 = sample2(index1);
	training2 = [training2part1 training2part2];
	id2 = [sample1ID(index1) sample2ID(index1)];
   
    for i = 1:size(id2,1)
    	pairBelong(id2(i,1), id2(i,2)) = 2;
    end

	%positive pair that belongs to group 3
	index1 = find(index == 3);
	training3part1 = sample1(index1);
	training3part2 = sample2(index1);
	training3 = [training3part1 training3part2];
	id3 = [sample1ID(index1)  sample2ID(index2)];

    %for pair 3 
	for i = 1:size(id3,1)
		pairBelong(id3(i,1), id3(i,2)) = 3;
	end

	%calculate the distance between negative pairs
	firstColumn = ux(:,negativePair(:,1));
	secondColumn = ux(:,negativePair(:,2));

	negPairDist1 = sum((Lmetric1 * (firstColumn - secondColumn)).^2);
	negPairDist2 = sum((Lmetric2 * (firstColumn - secondColumn)).^2);
	negPairDist3 = sum((Lmetric3 * (firstColumn - secondColumn)).^2);

	negPairDist = [negPairDist1 negPairDist2 negPairDist3];
	
	[largest, index] = max(negPairDist, [], 2);

	index1 = find(index == 1);
    for i = 1:size(index1,1)
		pairBelong(negativePair(index1(i),1), negativePair(index1(i),2)) = 1;
	end
	negativeGroup1 = [ux(:,negativePair(index1,1)); ux(:, negativePair(index1,2))];

	index1 = find(index == 2);
	for i = 1:size(index1,1)
		pairBelong(negativePair(index(i),1), negativePair(index2(i),2)) = 2;
	end
	negativeGroup2 = [ux(:,negativePair(index1,1)); ux(:, negativePair(index1,2))];

	index1 = find(index == 3);
	for i = 1:size(index1,1)
		pairBelong(negativePair(index(i),1), negativePair(index2(i),2)) = 3;
	end
	negativeGroup3 = [ux(:,negativePair(index1,1)); ux(:, negativePair(index1,2))];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




































