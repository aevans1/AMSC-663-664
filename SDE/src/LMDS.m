%Implementation of Landmark Multi-Dimensional Scaling for a given set of
%Landmarks and Data
function [embed_L,embed_Z] = LMDS(L,Z,rho,d)
%TODO: comments
	D = size(L,1); %each column of L is a landmark vector in R^D
	m = size(L,2); %number of landmark vectors
	p = size(Z,2); %number of data vectors
	
	%%%%Step 0:compute distance matrix for L%
	square_dist = zeros(m,m);
	for i = 1:m
		for j = 1:m
			square_dist(i,j) = ( rho(L(:,i),L(:,j)) )^2;	
		end
	end
	
	%%%Step 1: MDS for landmarks L
	%TODO: comment about output of this
	[embed_L,V,eigvals] = MDS(square_dist,rho,d);

	%Pseudo-inverse transpose of embed_L(verify? Taking from LMDS paper)
	pinv_L = (diag(eigvals))^(-1/2)*V.';
	
	%%%%Step 2: Project each data point to R^d based on Landmark embedding
	dist_z = ones(m,1); %initialize distance vectors for data
	embed_Z = zeros(d,p);  %initialize output list of embed data
	mean_L = mean(L,2); % mean landmark vector
	
	for i = 1:p
	%compute distance vector for z against landmarks, project onto embedded
	%...Landmark vectors
		for j = 1:m
			dist_z(j,1) = rho(Z(:,i),L(:,j)); 
		end
		embed_Z(:,i) = (-1/2)*pinv_L*(dist_z - mean_L);
	end

end

function [scaling,V,eigvals] = MDS(square_dist,d)
%perform Mult-Dimensional Scaling on data_set x_1,...x_m, with distance function rho, given a matrix of squared
%distances of the data points
%inputs: square-dist - m x m matrix of squared-distances of vectors x_1,..x_m
%					   i,jth entry of square dist is rho(x_i,x_j)^2
%		           d - desired dimension for embedding data set
%outputs: scaling - Matrix of embedded data points in R^d, columns are embedded
%					data vectors
%		        V - Matrix of first d eigenvectors of centered square-dist
%		            matrix, columns are eigenvectors
%		  eigvals - First d eigenvalues of centered square-dist matrix
	
	D = size(square_dist,1); %each column is length of data vector in R^D
	m = size(square_dist,2); %number of data vectors


	%%%%Center square_dist matrix
	H = (eye(m) - (1/m)*ones(m,1)*ones(m,1).'); %mean-centering matrix 
	B = (-1/2)*H*square_dist*H;

	%%%%Find eigvecs corresponding to highest d eigvals of B, sort
	%V is matrix of eigvecs of B
  	[V,eigvals] = eig(B);
  	[eigvals,ind] = sort(diag(eigvals),'descend');
	
	%only keep top d eigvals,corresponding eigvecs
	eigvals = eigvals(1:d);
  	V = V(:,ind(1:d));
	
	scaling =  diag(eigvals)^(1/2)*V.';

end
