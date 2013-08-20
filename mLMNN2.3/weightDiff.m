load('weightTest.mat');
y = {};
for j = 1:2
    uxLdim = selectSet{j}.metric;
    uxTest1 = test{j}.test1;
    uxTest2 = test{j}.test2;
    x = 1:316;
    y2 = [];
    for i = 1 : size(uxTest1,2)
        temp = bsxfun(@minus, uxLdim * uxTest2, uxLdim * uxTest1(:,i));
        temp = sum(temp.^2,1);
        
        temp1 = bsxfun(@minus, uxLdim * uxTest1, uxLdim * uxTest2(:,i));
        temp1 = sum(temp1.^2,1);

        temp = temp + temp1;
        [minimum, minIndex] = sort(temp);
        %myrank = myrank + minIndex;
        minIndex = find(minIndex == i);
        %y2 = y2 + [zeros(1,minIndex-1) ones(1,317-minIndex)];
        y2 = [y2 minIndex];
    end

    %Z = trapz(x,y2);
    %fprintf('%f\n',Z/(315*100));
    %plot(x,y2,'color','g');
    %hold on;
    %display(find(y2 == 1));
    y{j} = find(y2 == 1);
end