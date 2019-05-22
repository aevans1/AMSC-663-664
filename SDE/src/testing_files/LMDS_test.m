function LMDS_test()
%Verify that LMDS output is the projection of the data onto principal
%components of the landmarks
d = 2;
D = 7;

X = rand(D,10);
L = rand(D,3);
rho = @(x,y) norm(x - y);

%Center data landmarks
L = L - mean(L,2);
X = X - mean(X,2);

%%%Compute PCA of Landmarks
[proj_L,eigvecs_L,eigvals_L] = PCA(L,d);

%%%Compute LMDS of data X with landmarks L
[embed_L,embed_X] = LMDS(L,X,rho,d,false);
%%%Project data onto (up to d)eigenspace of L
PCA_proj = eigvecs_L'*X;
%%%Check that projection is same as LMDS
try
%Algorithms can find eigenvalues of opposite signs, need to allow for this
%in equivalence, using abs
    error = norm(abs(PCA_proj) - abs(embed_X))
    assert(error < 1e-4);
catch
    fprintf("LMDS not equivalent to PCA projection from landmarks! Check code \n");
end
end

function [projection,eigvecs,eigvals] = PCA(X,d)
%Input: columns of Y are feature vectors, epsilon is cutoff for principal
%component eigenvalues, lower than epsilon is not used
%output: D x N matrix, columns are data vectors
m = size(X,2);

Y = X - mean(X,2);

[eigvecs,eigvals] = eig((1/(m))*Y*Y.');

%sort eigenvectors in decreasing order of eigvals, keep
[eigvals,ind] = sort(diag(eigvals),'descend');
reduced_ind = ind(1:d);

%transform data
eigvecs = eigvecs(:,reduced_ind);
projection = (eigvecs.')*Y;
end



