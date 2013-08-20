%This program combines LBP and HSVYUV together as the features.
load('finalFeature.mat');
load('data/viper_features.mat');
mixFeature = [feature;ux(1:100,:)];

save('mixFeature.mat','mixFeature');