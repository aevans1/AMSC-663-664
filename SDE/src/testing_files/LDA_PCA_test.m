d =3;
D = 1;
m = 10;

X = randn(D,m);
rho = @(x,y) norm(x-y);
%%
%Y = X - mean(X,2);
Y = X;

%first, test MDS
%%%%Compute distance matrix for X%

%Z = Y - mean(Y,2);
square_dist = zeros(m,m);
for i = 1:m
	for j = 1:m
		square_dist(i,j) = ( rho(Y(:,i),Y(:,j)) )^2;	
	end
end

[mds_X,V,mds_eigvals] = MDS(square_dist,d);
[pca_X,eigvecs,pca_eigvals] = PCA(X,d);

% fprintf("Is pca(X) == mds_X?");
% pca_X == mds_X
 fprintf("pcaX ./ mdsx = ");
 pca_X./mds_X
% fprintf("pcaX - mdsX = ");
% 
% pca_X - mds_X
%mean(pca_X,1)
%mean(mds_X,2)
%mean(pca_X,2)
%mds_X - mean(mds_X,2)
%mds_X
%fprintf("original X: \n");

%fprintf("X ratio with PCA: \n");
%pca_X./X
%fprintf("X ratio with MDS: \n");
%mds_X./X
%%%fprintf("X difference with PCA: \n");
%%pca_X - X
%%fprintf("X difference with MDS: \n");
%%mds_X - X
%(mds_X + mean(Y,2))./Y
%pca_X./Y
%%mean(Y,2) + mds_X


%
%%
%%
%%
%%
%%
%L = X(:,1:2);
%L = L - mean(L,2);
%%%
%%%
%%%
%%%
%[PCA_L,eigvecs] = PCA(L,d);
%%%
%%test_X = (eigvecs.')*X;
%%%
%[new_L,new_X] = LMDS(L,X,rho,d);
%new_L - PCA_L
%new_L./L
%

%new_X - test_X
%
%
%

