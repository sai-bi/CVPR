%This program is used to do PCA
%setpaths
load('feature.mat');


outdim = 320;

HSV = [];
for i = 1:size(featureAHSV,2)
	HSV = [HSV featureAHSV{i}];
end

for i = 1:size(featureBHSV,2)
	HSV = [HSV featureBHSV{i}];
end
% 
 HSV(find(isnan(HSV))) = 0;
% hsvL = pca(HSV)';
% hsvL = hsvL(1:outdim,:);
% 
% 
YUV = [];
for i = 1:size(featureAYUV,2)
	YUV = [YUV featureAYUV{i}];
end

for i = 1:size(featureBYUV,2)
	YUV = [YUV featureBYUV{i}];
end

YUV(find(isnan(YUV))) = 0;
% 
% yuvL = pca(YUV)';
% yuvL = yuvL(1:outdim,:);
save('featureData','HSV','YUV');
%save('pcaData.mat','hsvL','yuvL');


hsvL = pca(HSV)';
hsvL = hsvL(1:outdim,:);
yuvL = pca(YUV)';
yuvL = yuvL(1:outdim,:);

HSVpart = hsvL * HSV;
YUVpart = yuvL * YUV;

feature = [HSVpart; YUVpart];

save('finalFeature.mat','feature','HSVpart','YUVpart');




















