load('dataLDA.mat');
X = [sample1 sample2];
L = pca(X)';
L = L(1:60,:);
sample1 = L * sample1;
sample2 = L * sample2;
[uxLdim,Det]= lmnn2([sample1 sample2],[1:316 1:316],1,'maxiter',15000,'checkup',0);
save('uxLdim2','uxLdim');
