load('multiExpert.mat');


%for j = 1:size(selectSet,2)
for j = 1:30
    uxLdim = selectSet{j}.metric;
    uxTest1 = selectSet{j}.test1;
    uxTest2 = selectSet{j}.test2;
    x = 1:316;
    y2 = zeros(1,316);
   % for i = 1 : size(uxTest1,2)
    for i = 1 : size(uxTest1,2)
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
    y2 = y2 /316 * 100;
    Z = trapz(x,y2);
    if(Z/(315*100) > 0.94)
        fprintf('%d: %f\n',j,Z/(315*100));
        plot(x,y2,'color','g');
        hold on;
    end
end