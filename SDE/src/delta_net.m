function [net,neighbors] = delta_net(init,delta,rho,is_random)
%sub-sample delta-net from given initial set of points, following the
%brute-force method given in section 3.1 of ATLAS paper
%inputs: init - D X N matrix, set of N vectors in R^d%
%		 delta - coarseness of delta-net,
%		 rho -  given distance function
%output: net - columns are data points, all distances >= delta apart
%		 neighbors - struct array, each net point index has a struct with key value
%		 'nbr' that stores the indices of close by net points
%NOTE: all points in domain should be within delta of any net point, not sure
%how to verify?
%
%
%
	%%%Default mode is random
	if ~exist('is_random','var'), is_random = true; end

	%%%%%%%%%%
	%Step 1: Build the Delta Net
	%%%%%%%%%%
	
	%%Random mode shuffles initial point set
	if is_random
		init = init(:,randperm(length(init)));
	end

	net = [];

	%%%check each point of init set, add to net if far enough away
	for n = 1:size(init,2)	
		if far(init(:,n),net,delta,rho)
			net(:,end+1) = init(:,n);
		end
		
		%fprintf("neighbors %d of %d initial points added \n",n,size(init,2));
	
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
		for m = 1:N
			if rho(net(:,n),net(:,m))<2*delta && m ~= n
				list(end+1) = m;
			end
		end
		neighbors(n).nbr = list;
		
		%fprintf("neighbors for net point %d of %d finished \n",n,N);
	
	end
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


