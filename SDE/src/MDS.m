function [scaling,V,eigvals] = MDS(square_dist,d)
%perform Mult-Dimensional Scaling on data_set X with distance function rho
%inputs: square-dist - m x m matrix of square-d distances from a data set X,
%with D rows and m columns, where each column is a data point
%	     D - input dimension of data vectors of X
%	     d - desired dimension for embedding data set
%outputs: scaling - Matrix of embedded data points in R^d, columns are embedded
%					data vectors
%		        V - Matrix of first d eigenvectors of centered square-dist
%		            matrix, columns are eigenvectors
%		  eigvals - First d eigenvalues of centered square-dist matrix
%NOTE: if rho is euclidean distance, this algorithm is exactly PCA on X	
	
	m = size(square_dist,1); %number of data vectors

	%%%%Center square_dist matrix
	H = (eye(m) - (1/m)*ones(m,1)*ones(m,1).'); %mean-centering matrix 
	B = (-1/2)*H*square_dist*H;

	%%%%Find eigvecs corresponding to highest d eigvals of B, sort
	%V is matrix of eigvecs of B
  	[V,eigvals] = eig(B);
            
  	[eigvals,ind] = sort(diag(eigvals),'descend');
	reduced_ind = ind(1:d);
	    
    %only keep top d eigvals,corresponding eigvecs
	eigvals = eigvals(1:d);
  	V = V(:,reduced_ind);
	scaling =  diag(eigvals)^(1/2)*V.';
end
