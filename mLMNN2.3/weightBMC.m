load('weightTest1.mat');

weightNum = size(selectSet,2);
weight = zeros(1,weightNum);
weightSum = 0;
z = -Inf;
randNum = 100;
likelihood = zeros(1,randNum);
dataSize = size(test{1}.test1,2);
for k = 1:randNum
	v = -log(unifrnd(0,1,1,weightNum));
	v = v/sum(v);

	x = 1:dataSize;
	y2 = zeros(1,dataSize);
	for i = 1 : dataSize
		%temp = zeros(1,dataSize);
	    myrank = zeros(1,dataSize);
	    for j = 1:size(selectSet,2)
	    	uxLdim = selectSet{j}.metric;
	    	uxTest1 = test{j}.test1;
	    	uxTest2 = test{j}.test2;
	    	temp2 = bsxfun(@minus, uxLdim * uxTest2, uxLdim * uxTest1(:,i));
			temp2 = sum(temp2.^2,1);
		    
		    temp1 = bsxfun(@minus, uxLdim * uxTest1, uxLdim * uxTest2(:,i));
			temp1 = sum(temp1.^2,1);

			temp =temp1 + temp2;   
			[minimum,minIndex] = sort(temp);
			for p = 1:dataSize
                myrank(p) = myrank(p) + find(minIndex == p) * v(j);
        	end      
	    end
	    [minimum, minIndex] = sort(myrank);
	    minIndex = find(minIndex == i);
		y2 = y2 + [zeros(1,minIndex-1) ones(1,dataSize + 1-minIndex)];
	end
	accu = y2(10)/dataSize;
	likelihood(k) = dataSize*(accu * log(accu) + (1-accu) * log(1-accu));

	if(likelihood(k) > z)
		weight = weight * (exp(z - likelihood(k)));
		z = likelihood(k);
	end
	w = exp(likelihood(k) - z);
	weight = weight * weightSum / (weightSum + w) + w * v;
	weightSum = weightSum + w;
    fprintf('Iter:%d, accu: %f\n',k,accu);
end

weight = weight/sum(weight);

save('weight.mat','weight');




















