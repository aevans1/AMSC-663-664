%Implemesntation of Landmark Multi-Dimensional Scaling for a given set of
%Landmarks and Data
function [embed_L,embed_Z] = LMDS(L,Z,rho,d,normalize)
%inputs:
%	L - set of landmarks for Z, array with columns in R^D
%		as landmark vectors
%	Z - set of data points, array with columns in $\R^D$
%		as data vectors
%	rho - given distance function
%	d - intrinsic dimension, LMDS projects columns of
%					Z,L onto R^d
%   normalize - boolean, true -> algorithm applies PCA to normalize output
%outputs:
%   embed_L - MDS output for landmarks, array with
%					columns in $\R^d$ as projected landmark vectors
%   embed_Z - projected $Z$ data, array with columns
%					in $\R^d$ as projected data vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%Adjust for extra parameters

%By default, normalize output via PCA
if ~exist('normalize','var')
    normalize = true; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%D = size(L,1); %each column of L is a landmark vector in R^D
	m = size(L,2); %number of landmark vectors
	p = size(Z,2); %number of data vectors

    
    meanvec = mean(L,2);
    L = L - meanvec;
    Z = Z - meanvec;
        
	square_dist = compute_square_dist(L,rho);
   
	%%%Step 1: MDS for landmarks L
	%TODO: comment about output of this
	[embed_L,V,eigvals] = MDS(square_dist,d);
        
    %Pseudo-inverse transpose of embed_L
	pinv_L = (diag(eigvals))^(-1/2)*V.';
	%%%%Step 2: Project each data point to R^d based on Landmark embedding
	dist_z = ones(m,1); %initialize distance vectors for data
	embed_Z = zeros(d,p);  %initialize output list of embed data
	mean_square_dist = mean(square_dist,2); % mean landmark vector

	for i = 1:p
	%compute distance vector for z against landmarks, project onto embedded
	%Landmark vectors
		for j = 1:m
			dist_z(j,1) = rho(Z(:,i),L(:,j))^2; 
		end
		embed_Z(:,i) = (-1/2)*pinv_L*(dist_z - mean_square_dist);
    end
    
    
    if normalize
        %%%Step 3: Normalize data by applying PCA
        data = [embed_L embed_Z];
        PCA_data = PCA(data,d);
        embed_L = PCA_data(:,1:size(L,2));
        embed_Z = PCA_data(:,size(L,2) + 1 : end);
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
	    
    %only keep top d positive eigvals,corresponding eigvecs
	eigvals = eigvals(1:d);
    eigvals = eigvals(eigvals > 0);
    if isempty(eigvals)
        fprintf("landmark sq. dist matrix has no positive eigvals!");
    end
%     [lam,isort] = sort(diag(E),'descend');
% V = V(:,isort);
% ind = find(lam > 0);
% if isempty(ind)
%     fprintf('The landmark distance squared matrix has no positive eigenvalues\n');
%     return;
% end
% kplus = min(k,length(ind));
% lamk = lam(1 : kplus);

    
    
    
  	V = V(:,reduced_ind);
    size(eigvals)
    size(V.')
	scaling =  diag(eigvals)^(1/2)*V.';
end

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
