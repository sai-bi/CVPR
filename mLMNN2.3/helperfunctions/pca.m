function [evects,evals] = pca(X)
% [evects,evals] = pca(X)
%e
% finds principal components of 
%
% input: 
%  X  dxn matrix (each column is a dx1 input vector)
% 
% output: 
% evects  columns are principal components (leading from left->right)
% evals   corresponding eigenvalues
%
% See also applypca
%
% copyright by Kilian Q. Weinberger, 2006

%*************************************************************************

% d: number of dimensions
% N: number of samples
[d,N]  = size(X);

X=bsxfun(@minus,X,mean(X,2));  
cc = cov(X',1); % compute covariance 
save('cc.mat','cc');
[cvv,cdd] = eig(cc); % compute eignvectors
[zz,ii] = sort(diag(-cdd)); % sort according to eigenvalues
evects = cvv(:,ii); % pick leading eigenvec
cdd = diag(cdd); % extract eigenvalues
evals = cdd(ii); % sort eigenvalues


