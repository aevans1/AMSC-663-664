function [square_dist] = compute_square_dist(X,rho)
%X a square matrix, computes squared distances between all columns

%%%%Compute distance matrix for X%
	square_dist = zeros(size(X,2),size(X,2));
	for i = 1:size(X,2)
		for j = 1:size(X,2)
			square_dist(i,j) = ( rho(X(:,i),X(:,j)) )^2;	
		end
    end
end
