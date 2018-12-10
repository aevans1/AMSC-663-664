%Input: columns of Y are feature vectors, epsilon is cutoff for principal
%component eigenvalues, lower than epsilon is not used
%inputs: D x N matrix, columns are data vectors
function [projection,eigvecs,eigvals] = PCA(X,d)


Y = X - mean(X,2);
[eigvecs,eigvals] = eig(cov(Y.',1));

%sort eigenvectors in decreasing order of eigvals, keep
[eigvals,ind] = sort(diag(eigvals),'descend');
reduced_ind = ind(1:d);
%eng = (1/sum(eigvals))*cumsum(eigvals.');
%reduced_ind = (eng < energy).*ind.';
%reduced_ind = reduced_ind(reduced_ind > 0);
%dim = length(reduced_ind);

%transform data
eigvecs = eigvecs(:,reduced_ind);
projection = (eigvecs.')*X;
end
