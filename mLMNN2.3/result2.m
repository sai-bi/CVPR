load('dataLDA.mat');
load('active.mat');
load('uxLdim2.mat');
uxTest1 = uxTest1' * result;
uxTest2 = uxTest2' * result;

uxTest1 = uxTest1';
uxTest2 = uxTest2';

x = 1:316;
y2 = zeros(1,316);


for i = 1:316
    temp = bsxfun(@minus, uxLdim * uxTest2, uxLdim * uxTest1(:,i));
    temp = sum(temp.^2,1);

    temp1 = bsxfun(@minus, uxLdim * uxTest1, uxLdim * uxTest2(:,i));
    temp1 = sum(temp1.^2,1);
    
    temp = temp + temp1;
    [minimum, minIndex] = sort(temp);
    %myrank = myrank + minIndex;
    minIndex = find(minIndex == i);
    y2 = y2 + [zeros(1,minIndex-1) ones(1,317-minIndex)];
end
y2 = y2*100/316;
plot(x,y2,'color','g');
hold on;