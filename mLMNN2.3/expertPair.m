load('expertData.mat');
load('expertLdim.mat');

metric1 = Ldim;
metric2 = Ldim;
metric3 = Ldim;

[x,y] = meshgrid(1:316, 317:632);
negativePair = [x(:) y(:)];         %generate all pairs
negativePair(1:316:(316*316),:) = []; %remove positive pair

pairPerm = randperm(99540); 
pairDist = zeros(632,632); %used to record the distance between each pair
%random partition

for i = 1:33180
	point1 = sample(:,negativePair(pairPerm(i),1));
	point2 = sample(:,negativePair(pairPerm(i),2));
	pairDist(negativePair(pairPerm(i),1), negativePair(pairPerm(i),2)) = sum((metric1*(point1 - point2)).^2);
end

for i = 33181:66360
	point1 = sample(:,negativePair(pairPerm(i),1));
	point2 = sample(:,negativePair(pairPerm(i),2));
	pairDist(negativePair(pairPerm(i),1), negativePair(pairPerm(i),2)) = sum((metric2*(point1 - point2)).^2);
end

for i = 66361:99540
	point1 = sample(:,negativePair(pairPerm(i),1));
	point2 = sample(:,negativePair(pairPerm(i),2));
	pairDist(negativePair(pairPerm(i),1), negativePair(pairPerm(i),2)) = sum((metric3*(point1 - point2)).^2);
end

save('negPair.mat','pairDist','negativePair');