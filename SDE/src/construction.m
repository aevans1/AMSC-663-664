function [new_Sim,neighbors,net] = construction(S,init,delta,rho,m,p,t_0,d,load_net)
%Constructs delta-net, landmarks, ATLAS, and SDE simulator
%inputs: S - SDE simulator
%        init - D X N matrix, set of N vectors in R^d
%        delta - homogenization scale, affects density of sample net
%		 rho -  given distance function
%		 m - desired number of landmarks to generate for each net point
%        p - num sample paths for each point in net
%		 t_0 - desired time of simulation for each call to S
%        d - intrinsic dimension of dynamics, problem specifics
%        load_net - boolean, true value loads up atlas variables specified in
%        			'preload.m', such as a delta net and neighbors
%outputs: new_Sim - struct with key values that specify parameters for
%                   simulating a constant coefficient SDE 
%         net - columns are data points, all distances >= delta apart
%  		  neighbors - struct array, each net point index has a struct with key value
%		 'nbr' that stores the indices of close by net points

	
%%%Create delta_net from initial points, create landmarks for delta_net
%	TODO: better way to do this? Loading up structs if saved
	if ~exist('load_net','var')
		load_net = false;
	end
	
	if load_net == false  
		[net,neighbors] = delta_net(init,delta,rho);
	else
		fprintf("Using loaded variables for delta net and neighbors\n");
		load('preload','net','neighbors');
		fprintf('net is...\n');
		net
	end

	%fprintf("delta net is.... \n");
	%net
	
	D = size(net,1); %assuming each column of delta_net is a data_point
	N = size(net,2); %number of points in delta_net
	
	%%%Initialize arrays, and create landmarks from the net		
	X = zeros(D,p,N);
	L(N).L = {};
	embed_L(N).L = {}; %LMDS embeddings of landmarks
	embed_X(N).X = {}; %LMDS embeddings of simulating points
	Phi(N).Phi =[]; 	%multipliers for 1-dim chart mappings
	
	T(N,N).T = {}; %T is the collection of transition maps from chart to chart
	C(N,N).C = {}; %C is the collection of chart centers for each chart
	mu(N,N).mu = {}; %mu is the collection of average landmark for each chart

	%%%Create landmarks
	A = zeros(D,m+1,N); %landmark set
	for n = 1:N
	%generate m paths, save endpoints for delta_net point n
	%NOTE: %first point in landmark list is the %delta_net point itself
			A(:,1,n) = net(:,n);
			A(:,2:m+1,n) = S(net(:,n),m,t_0);
	end
	
	%%%Create atlas
	for n = 1:N
		fprintf("Constructing for net point %d \n",n);
		
		%simulate p paths around net point y_n
		X(:,1:p,n) = S(net(:,n),p,t_0);

		%%%%Take union of all neighboring landmarks to y_n
		%net_nbr contains global indices of neighbor net points to netpoint n	
		net_nbr = neighbors(n).nbr;
		num_nbr = length(net_nbr);

		%store points from neighboring landmarks as collection of columns
		%associated to y_n
		L(n).L = reshape(A(:,:,net_nbr),D,(m+1)*num_nbr);
		
		%%%Step 1:Create chart for each net point
		[local_L,local_X,local_Phi] = create_chart(X(:,:,n),L(n).L,rho,d,D);
		
		embed_L(n).L = local_L;
		embed_X(n).X = local_X;
		Phi(n).Phi = local_Phi;

		%%%Step 2: constructing local SDE
		[T,C,B(:,n),Sigma(:,:,n),mu] = construct_local_SDE(n,neighbors,embed_L,embed_X,p,t_0,m,d,T,C,mu);

	end

	new_Sim.T = T;
	new_Sim.B = B;
	new_Sim.C = C;
	new_Sim.Sigma = Sigma;
	new_Sim.mu = mu;
	
	%matters for D = d = 1 case only!
	new_Sim.Phi = Phi;	
	
	filename =[datestr(now, 'dd_mmm_yyyy_HH_MM'),'_','d',num2str(d),'_','construction'];

	save(filename);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [local_L,local_X,Phi] = create_chart(X,L,rho,d,D);
%TODO:comments	
      			
		%tic
		[local_L,local_X] = LMDS(L,X,rho,d);
        %LMDS_time = toc

		Phi = [];
		if d ==1 && D ==1
		%%%%Note:D = d = 1 case only!
		%Here Phi is whatever was multiplied by X's to get embedding
			temp_X = X - mean(L,2); 
			Phi = local_X(:,1)/temp_X(:,1);
		end
		%%%%%

        %%%Centering all data for chart around chart center
		center = local_L(:,1);
		local_L = local_L - center;
		local_X = local_X - center;
	end

function [T,C,B,Sigma,mu] = construct_local_SDE(n,neighbors,embed_L,embed_X,p,t_0,m,d,T,C,mu)
%TODO: comments

	%net_nbr contains global indices of neighbor net points to netpoint n	
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
	B = (1/ (p*t_0) )*sum(local_X,2); %drift for y_n

	%fprintf("local average before normalizing for time: %f \n",(1/p)*sum(local_X,2));
	Sigma = sqrt(1/t_0)*sqrtm(cov(local_X.')); %diffusion for y_n


	%%%Compute switching maps
		
	%NOTE: net_nbr = [n idx1 idx 2 ...], idx1 < idx 2 < ...
	for i = 1:num_nbr	
		
		if net_nbr(i) < n
			j = net_nbr(i);	%global index of neighbor of n

			%L_j' in atlas paper
			local_nbr_L = embed_L(j).L; %embeddings wrt nbr of y_n

			%NOTE: indexing to find nbr is based on that net_nbr = [n nbr1 nbr2
			%.] [y_n landmarks| y_nbr1 landmarks ...]

			%From paper: this is L_k,j
			%embedded landmarks pair wrt y_n
			L_nj(:,1:m+1) = local_L(:,1:m+1);  %y_k wrt y_k embedding	
			L_nj(:,m+2:2*m+2) = local_L(:, (i - 1)*(m+1) + 1: i*(m+1)); %nbr wrt y_k embedding
			
			mu(n,j).mu = mean(L_nj,2);

			%%%find y_n wrt nbr embedding
			nbr_switch_ind = neighbors(j).nbr;
			ind = find(nbr_switch_ind ==n);

			%From paper: this is L_j,k
			L_jn(:,1:m+1) = local_nbr_L(:,1:m+1);  %nbr i wrt nbr i embedding	
			L_jn(:,m+2:2*m+2) = local_nbr_L(:,(ind-1)*(m+1) + 1: ind*(m+1)); %y_n  wrt nbr i embedding

			mu(j,n).mu = mean(L_jn,2);
			
			%T(n,j).T is transition map from chart n to chart j
			%i.e, T(n,j).T*x changes x from chart n coords to chart j
			%coords
 			T(n,j).T = (L_jn - mu(j,n).mu)*pinv(L_nj - mu(n,j).mu);
			T(j,n).T = (L_nj - mu(n,j).mu)*pinv(L_jn - mu(j,n).mu);
		end
		%end transition map computing for neighbor
    end
		%end computation of transition maps
end
