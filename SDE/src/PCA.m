%Input: columns of Y are feature vectors, epsilon is cutoff for principal
%component eigenvalues, lower than epsilon is not used
%inputs: D x N matrix, columns are data vectors
function [projection,eigvecs,eigvals] = PCA(X,d)

m = size(X,2);

Y = X - mean(X,2);

[eigvecs,eigvals] = eig((1/m)*Y*Y.');

%sort eigenvectors in decreasing order of eigvals, keep
[eigvals,ind] = sort(diag(eigvals),'descend')
reduced_ind = ind(1:d);

%transform data
eigvecs = eigvecs(:,reduced_ind);
projection = (eigvecs.')*X;

end
