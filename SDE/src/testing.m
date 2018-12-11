d =1
m = 4

X = randn(d,m);
rho = @(x,y) norm(x-y);

X = X - mean(X,2);
% = X - mean(X,1);

%first, test MDS
%%%%Compute distance matrix for X%
square_dist = zeros(m,m);
for i = 1:m
	for j = 1:m
		square_dist(i,j) = ( rho(X(:,i),X(:,j)) )^2;	
	end
end

%%%%Center square_dist matrix
H = (eye(m) - (1/m)*ones(m,1)*ones(m,1).'); %mean-centering matrix 
B = (-1/2)*H*square_dist*H



[mds_X,V,pca_eigvals] = MDS(square_dist,d);

[pca_X,eigvecs,mds_eigvals] = PCA(X,d);

%(X.'*X) - B

fprintf("size mds: %d \n",size(mds_X));
size(pca_X)

%size(pca_X)
pca_X./X
mds_X./X

pca_eigvals
mds_eigvals

%
%
%
%
%
L = X(:,1:2);
L = L - mean(L,2);
%%
%%
%%
%%
[PCA_L,eigvecs] = PCA(L,d);
%%
%test_X = (eigvecs.')*X;
%%
[new_L,new_X] = LMDS(L,X,rho,d);
new_L - PCA_L
new_L./L


%new_X - test_X
%
%
%

