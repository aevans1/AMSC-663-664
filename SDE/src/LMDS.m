%Implementation of Landmark Multi-Dimensional Scaling for a given set of
%Landmarks and Data
function [embed_L,embed_Z] = LMDS(L,Z,rho,d)
%TODO: comments
	%D = size(L,1); %each column of L is a landmark vector in R^D
	m = size(L,2); %number of landmark vectors
	p = size(Z,2); %number of data vectors

	%TESTING:
	%L = L - mean(L,2);
	%Z = Z - mean(Z,2);

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
	%...Landmark vectors
		for j = 1:m
			dist_z(j,1) = rho(Z(:,i),L(:,j))^2; 
		end
		embed_Z(:,i) = (-1/2)*pinv_L*(dist_z - mean_square_dist);
	end
 	embed_Z = PCA(embed_Z,d);
end


