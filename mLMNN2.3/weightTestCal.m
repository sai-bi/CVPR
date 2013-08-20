load('weightTest2.mat');
load('weight.mat');

%weight = ones(1,size(selectSet,2));
dataSize = size(test{1}.test1,2);
x = 1:dataSize;
y2 = zeros(1,dataSize);
for i = 1 : dataSize
	%temp = zeros(1,dataSize);
    myrank = zeros(1,dataSize);
    for j = 1:size(selectSet,2)
        temp = zeros(1,dataSize);
    	uxLdim = selectSet{j}.metric;
    	uxTest1 = test{j}.test1;
    	uxTest2 = test{j}.test2;
    	temp2 = bsxfun(@minus, uxLdim * uxTest2, uxLdim * uxTest1(:,i));
		temp2 = sum(temp2.^2,1);
	    
	    temp1 = bsxfun(@minus, uxLdim * uxTest1, uxLdim * uxTest2(:,i));
		temp1 = sum(temp1.^2,1);

        temp = temp + (temp1 + temp2);  
        [minimum, minIndex] = sort(temp);    
        for p = 1:dataSize
            myrank(p) = myrank(p) + find(minIndex == p) * weight(j);
        end  
    end
    [minimum, minIndex] = sort(myrank);

    minIndex = find(minIndex == i);
    %display(minIndex);
	y2 = y2 + [zeros(1,minIndex-1) ones(1,dataSize + 1-minIndex)];
end

y2 = y2/dataSize * 100;

Z = trapz(x,y2);
fprintf('%f\n',Z/((dataSize - 1)*100));
plot(x,y2,'color','g');
hold on;
