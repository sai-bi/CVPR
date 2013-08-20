setpaths
load('expertData.mat');
[Ldim,Det]=multiLMNN(sample,[1:316 1:316],1,'quiet',1,'maxiter',100,'checkup',0,'outdim',34);
save('expertLdim','Ldim');