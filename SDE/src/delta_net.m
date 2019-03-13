function [net,neighbors] = delta_net(init,delta,rho,is_random)
%sub-samples initial set of points to create a delta-net, following the
%brute-force method given in section 3.1 of ATLAS paper
%inputs: init - D X N matrix, set of N vectors in R^d%
%		 delta - coarseness of delta-net,
%		 rho -  given distance function
%output: net - columns are data points, all distances >= delta apart
%		 neighbors - struct array, each net point index has a struct with key value
%		 'nbr' that stores the indices of close by net points
%
%
%
	%%%Default mode is random
	if ~exist('is_random','var'), is_random = true; end

	%%%%%%%%%%
	%Step 1: Build the Delta Net
	%%%%%%%%%%
	
	%%%Random mode shuffles initial point set
	if is_random
		init = init(:,randperm(length(init)));
	end

	net = [];

	%%%check each point of init set, add to net if far enough away
	n = 1;
	init_copy = init;
	while n <= size(init_copy,2)
		if far(init_copy(:,n),net,delta,rho)
			net(:,end+1) = init_copy(:,n);
			init_copy(:,n) = [];
		else
			n = n+1;
		end
		
		%fprintf("neighbors %d of %d initial points added \n",n,size(init,2));
	end

	%%%Check that there are no isolated net points
	%NOTE: this does not check for if the delta net adequately covers the
	%domain
	for n = 1:size(net,2)
		if ~close(net(:,n),net,delta,rho)
			fprintf("isolated point!consider re-initializing delta_net \n");
		end
	end

	%%%%%%%%%%
	%Step 2: Find and store neighbors of each delta net point
	%%%%%%%%%%

	%%%Create struct of neighbors for delta-net: two net points y_n,y_m are connected if
	%rho(y_n,y_m) < 2delta
	N = size(net,2);
	neighbors = {};
	for n = 1:N
		list =[n]; %first neighbor of y_n will be y_n itself

		%%%check net for neighbors, add those indices to struct
		num_nbr = 0;
		for m = 1:N
			if rho(net(:,n),net(:,m))<2*delta && m ~= n
				list = [list; m];
			end
		end
		neighbors(n).nbr = list;
		
		%fprintf("neighbors for net point %d of %d finished \n",n,N);
	
	end
		%TODO: figure out whether to add in this covering check or not	
		%if num_nbr == 1:
		%	if is_random:	
		%		%create new point in random direction delta to 2 delta away
		%		net(:,end+1) = net(:,n) + (delta*rand() + delta)*rand(D,1);
		%	else:
		%		%create new point by adding 1.5 delta[ 1 1 1 1 1 ...]^T
		%		net(:,end+1) = net(:,n) + delta;
		%end
		

end

%%NOTE: cross-ref with test function file
%%%NOTE: if input x is part of net, this will return false
function is_far = far(x,net,delta,rho)
%determine if point 'x' is at least delta-far from set of points 'net'
%inputs: x is vector in R^D, net is D X N matrix, N vectors in net
%output: boolean, true if point is far from set
	is_far = true;
	N = size(net,2);
	k = 1;

	%Check is x is far away from net
	while (is_far && k <= N)
		if rho(x,net(:,k)) < delta
			is_far = false; %break!
		else
			k = k+1;
		end
	end
end

function is_close = close(x,net,delta,rho)
%determine if point 'x' is closer than 2delta to some point in net
%%inputs: x is vector in R^D, net is D X N matrix, N vectors in net
%output: boolean, true if point is far from set
	is_close = false;
	N = size(net,2);
	k = 1;

	%Check is x is far away from net
	while (~is_close && k <= N)
		if rho(x,net(:,k)) > 0 & rho(x,net(:,k)) < 2*delta 
			is_close = true; %break!
		else
			k = k+1;
		end
	end
end


