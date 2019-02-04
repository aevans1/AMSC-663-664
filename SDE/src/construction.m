function [new_Sim,neighbors,net] = construction(S,init,delta,rho,m,p,t_0,d)
%Constructs delta-net, landmarks, ATLAS, and SDE simulator
%inputs: S - SDE simulator
%         init - D X N matrix, set of N vectors in R^d
%        delta - homogenization scale, affects density of sample net
%		 rho -  given distance function
%		 m - desired number of landmarks to generate for each net point
%        p - num sample paths for each point in net
%		 t_0 - desired time of simulation for each call to S
%        d - intrinsic dimension of dynamics, problem specifics
%outputs: new_Sim - struct with key values that specify parameters for
%                   simulating a constant coefficient SDE 
%         net - columns are data points, all distances >= delta apart
%  		  neighbors - struct array, each net point index has a struct with key value
%		 'nbr' that stores the indices of close by net points

	
%%%Create delta_net from initial points, create landmarks for delta_net

	%TEMPORARY
	%fprintf("WARNING: temporarily using loaded variables for delta net and
	%neighbors \n");
	%load('atlas_driver','net','neighbors');
	
	[net,neighbors] = delta_net(init,delta,rho);
	
	fprintf("delta net is.... \n");
	net
	A = create_landmarks(S,net,m,t_0);
	
	D = size(net,1); %assuming each column of delta_net is a data_point
	N = size(net,2); %number of points in delta_net
	X = zeros(D,p,N);
	%nbr_landmarks = {}; %currently using a struct for collecting neighbor landmarks
	%NOTE: should join this struct with the neighbors struct!
		
	%NOTE: better way to write this?
	%For now, calculating the embedding for each data point, then doing ANOTHER
	%for loop for rest of algorithm
	for n = 1:N
		%simulate p paths around net point y_n
		X(:,1:p,n) = S(net(:,n),p,t_0);

		%NOTE: redundancy here, over collecting landmarks?
		%Take union of all neighboring landmarks to y_n
		net_nbr = neighbors(n).nbr;
		num_nbr = length(net_nbr);

		%store  points from neighboring landmarks as collection of columns
		%associated to y_n
		L = reshape(A(:,:,net_nbr),D,(m+1)*num_nbr);

        %Comparing mean-0 versions of L and X with LMDS outputs, 1-d only
        temp_L = L - mean(L,2);
		temp_X = X(:,:,n) - mean(X(:,:,n),2);		
       
				%TESTING: timing LMDS with tic/toc
		%tic
		%%Still don't know which of these to pick
		[embed_L(n).L,embed_X(n).X] = LMDS(L,X(:,:,n),rho,d);
	   %[embed_L(n).L,embed_X(n).X] = LMDS(temp_L,X(:,:,n),rho,d);
        %LMDS_time = toc

		local_L = embed_L(n).L;
		local_X = embed_X(n).X;
 		old_X = X(:,:,n);

		%TODO: change this
		Phi(n).Phi = [];
		%%%%Note:D = d = 1 case only!
		%Here Phi is whatever was multiplied by X's to get embedding
		if d ==1 && D ==1
			Phi(n).Phi = local_X(:,1)/temp_X(:,1);
		end
		%Phi(n).Phi = local_X(:,1)/old_X(:,1);
		%Phi(n).Phi = local_L(:,1)/temp_L(:,1);
		%%%%%

        %TESTING: various comparisons to see what LMDS is doing in 1-D
% 		fprintf("ratio with original L: \n")
% 		local_L(:,1)./L(:,1)
% 
%  		fprintf("ratio with mean-centered L: \n")
%  		local_L(:,1)./temp_L(:,1)		
% 		
   		%fprintf("ratio with original X: \n")
   		%local_X(:,1)/old_X(:,1)
% 
% 		fprintf("ratio with mean-centered X: \n")
% 		local_X(:,1)./temp_X(:,1)

		%TESTING: seeing what happens with no LMDS for 1d	
		%embed_L(n).L = L*Phi(n).Phi;
		%embed_X(n).X = X(:,:,n)*Phi(n).Phi;

        %TESTING: collecting non-centered chart images of net points
        %image(n).image = local_L(:,1);
        
		%%Centering all data for chart around chart center
		embed_L(n).L = embed_L(n).L - local_L(:,1);
		embed_X(n).X = embed_X(n).X - local_L(:,1);

	end

	for n = 1:N
	
    %TODO: write function to handle this for loop
    %construct_SDE(n,neighbors,C,B,Sigma,local_L,local_X,embed_L,embed_X,m,p,t_0)
	
    %%%%%Compute new simulator
		net_nbr = neighbors(n).nbr;
		num_nbr = length(net_nbr);
		
		local_L = embed_L(n).L;
		local_X = embed_X(n).X;

		%select neighboring delta-net points:
		%local_L = [y_n |y_n landmarks |nbr1 |nbr 1landmarks |nbr2|nbr2
		%landmarks...];
		ind = [1:m+1:(m+1)*num_nbr];

		%%%Compute chart 'C_n' for y_n and center around embedded y_n
		for i = 1:num_nbr
			j = net_nbr(i); %global index of neighbor
			C(n,j).C = local_L(:,ind(i)); %find nbr and center
		end
	
		%%%Compute Diffusion coefficients, drift coefficients around y_n
		%fprintf("diffusion coefficient for delta-net point %f \n",net(n));
		B(:,n) = (1/ (p*t_0) )*sum(local_X,2); %drift for y_n

		%fprintf("local average before normalizing for time: %f \n",(1/p)*sum(local_X,2));
		Sigma(:,n) = sqrt((1/t_0)*cov(local_X.',1)); %diffusion for y_n
		
		%%%Compute switching maps
		%Grab max index of nbr with global index less than n
		%Switch ind corresponds to indices of L matrix rather than delta_net
		
        switch_ind = net_nbr(net_nbr < n);
	
		%NOTE: net_nbr = [n idx1 idx 2 ...], idx1 < idx 2 < ...
		%NOTE: change this, too confusing, but switch_ind+1 is index of y_n
		%landmarks in L matrix
		for i = 1:length(switch_ind)
			
			j = switch_ind(i);
			local_nbr_L = embed_L(j).L; %embeddings wrt nbr of y_n

			%NOTE: indexing to find nbr is based on that net_nbr = [n idx1 idx2
			%.] and switch_ind = [idx1 idx2 ...], switch_ind + 1 will be index
			%of neighbor in landmakrs, [y_n landmarks| idx1 landmarks ...]
			ind = i + 1;

			%From paper: this is L_k,j
			%embedded landmarks pair wrt y_n
			L_nj(:,1:m+1) = local_L(:,1:m+1);  %y_k wrt y_k embedding	
			L_nj(:,m+2:2*m+2) = local_L(:,(ind-1)*(m+1):(ind)*(m+1) - 1); %nbr wrt y_k embedding

			mu(n,j).mu = mean(L_nj,2);

			%%%find y_n wrt nbr embedding
			%NOTE: re-write this
			nbr_switch_ind = neighbors(j).nbr;
			ind = find(nbr_switch_ind ==n);

			%From paper: this is L_j,k
			L_jn(:,1:m+1) = local_nbr_L(:,1:m+1);  %nbr i wrt nbr i embedding	
			L_jn(:,m+2:2*m+2) = local_nbr_L(:,(ind-1)*(m+1): ind*(m+1) - 1); %y_n  wrt nbr i embedding

			mu(j,n).mu = mean(L_jn,2);

			%T(n,j).T is transition map from chart n to chart j
			%i.e, T(n,j).T*x changes x from chart n coords to chart j
			%coords
            
 			T(n,j).T = (L_jn - mu(j,n).mu)*pinv(L_nj - mu(n,j).mu);
 			T(j,n).T = (L_nj - mu(n,j).mu)*pinv(L_jn - mu(j,n).mu);

        end
	end
	
	%%%Simulator now computed!	
	%TESTING
% 	for n = 1:N
% 		fprintf("neighbors of %d: \n",n);
% 		nbr = neighbors(n).nbr
% 		fprintf("transitions from index %d \n",n);
% 		for i = 1:length(nbr)
% 			fprintf("%d to %d \n",n,nbr(i));
% 			T(n,nbr(i)).T	
% 		end
% 	end
% 	
	%Now: we have chart centers C, drift coef B, diffusion coeff Sigma
	new_Sim.T = T;
	new_Sim.B = B;
	new_Sim.C = C;
	new_Sim.Sigma = Sigma;
	new_Sim.mu = mu;
	
	%matters for D = d = 1 case only!
	new_Sim.Phi = Phi;	
	%
    
    %TESTING
    %new_Sim.image = image;
end

