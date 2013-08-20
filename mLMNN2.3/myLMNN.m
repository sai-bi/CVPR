function [metric]= myLMNN(metric,pairDist,id,training,trainingLabel,maxiter, sample)
% Parameter meaning:
% metric: inital value for L
% metric1: distance metric for group 1
% metric2: distance metric for group 2
% metric3: distance metric for group 3
% pairDist: the distance of a negative pair
% id: the id of each data point
% training: training label, each column is a feature
% trainingLabel: a row vector, the label of each feature
% maxiter: maximum number of iteration
% This is a modified version of LMNN implementation.

[row, column] = size(training);
display(row);
display(column);
training1 = training(:,1:column/2);
training2 = training(:,(column/2+1) : end);

%calculate the outer product of positve pairs, that is, sum(Cij),i is a target neighbour
%of j
temp = training1 - training2;
sumCij = zeros(631, 631);

for i = 1 : size(training1,2)
	sumCij = sumCij + temp(:,i) * (temp(:,i)');
end

%some parameter settings
pars.stepsize = 1e-07;
pars.stepgrowth = 1.01;
pars.minstepsize = 0;
stepsize = pars.stepsize;

obj = zeros(1,maxiter);
objOld = [];
metricOld = zeros(34,631);

for iter  = 1 : maxiter
	%calculate the second part of Gt, that is Cijl, j is a target neighbour of i,
	%and l is an impostor of i and j.
    fprintf('Iteration %d=============================\n', iter);    
	sumCijl = zeros(631, 631);
	constraint = 0;
    
    %save old metric and obj value
    objOld = obj(max(iter-1,1));
    metricOld = metric;
    
    
    
	for i = 1:size(training1,2)
		pointID = id(i);
		pairDistance = pairDist(pointID,:); 

        posPair = metric * (training1(:,i) - training2(:,i));
        distIJ = posPair' * posPair;
		
		index = find(pairDistance < distIJ + 1);
		constraint = constraint + size(index,1);
        Cij = posPair * posPair';
        
		for l = 2 : size(index,1)
			negPair = training1(:,i) - sample(:,index(l)+316);
			Cijl = negPair * negPair';
			sumCijl = sumCijl + Cij - Cijl;
		end
    end

    
    
	Gt = sumCij * 0.5 + sumCijl * 0.5;
	
    if(iter > 1)
       metric=step(metric,Gt,stepsize,pars);
    end

	%calculate object function
	objPart1 = trace((0.5 * sumCij) * metric' * metric);
	objPart2 = trace((0.5 * sumCijl) * metric' * metric);
    objPart3 = constraint * 0.5;

    obj(iter) = objPart1 + objPart2 + objPart3;

    delta=obj(iter)-obj(max(iter-1,1));


    if(i>1 && delta>0) 
	 stepsize=stepsize*0.5;
     obj(iter) = objOld;
     metric = metricOld;
	 fprintf('***correcting stepsize***  \ndelta: %f\n', delta);
	 if(stepsize<pars.minstepsize) 
	 	stepsize=pars.minstepsize;
	 	fprintf('stepsize is 0 now :(');
	 end
	else 
	  stepsize=stepsize*pars.stepgrowth;
    end
    
   
end


fprintf('LMNN has finished ^_^ ');
save('metric.mat', 'metric');






function L=step(L,G,stepsize,pars)
% do step in gradient direction
G=2.*(L*G); 
L=L-stepsize.*G;     
return;
end

end
