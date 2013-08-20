load('uxData.mat');
load('uxLdim2.mat');

posPair = abs(sample1 - sample2);
posPair = abs(uxLdim * posPair);

negPair = zeros(631,99540);
column = 1;
for i = 1 : 316
    fprintf('%d ++++++++++\n',i);
	for j = 1 : 316
		if i == j
			continue;
		end
		temp = abs(sample1(:,i) - sample2(:,j));
		negPair(:,column) = temp;
        column = column+1;
	end
end

negPair = abs(uxLdim * negPair);

display(size(negPair));

save('posData.mat', 'posPair');
save('negData.mat', 'negPair');