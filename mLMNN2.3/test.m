load('negData.mat');
load('posData.mat');

split = zeros(1,60);
percent = zeros(1,316);
percent2 = zeros(1,316);
myrank = zeros(1,60);
for i = 1:1
	pos = posPair(i,:);
	neg = negPair(i,:);

	pos = sort(pos);
	neg = sort(neg);

	for j = 1:316
		index = find(pos < pos(1,j));
        index2 = find(neg < pos(1,j));
		percent(1,j) = size(index,2)/316.0;
        percent2(1,j) = size(index2,2)/95440.0;
    end
    x = linspace(1,316,316);
    plot(pos,percent,'r');
    hold on;
    plot(pos,percent2,'g');
    hold on;
end

