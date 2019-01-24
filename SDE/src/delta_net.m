function [net,neighbors] = delta_net(init,delta,rho)
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
	k = 1; %first net point index
	%k = randi(length(init)); %optional: random initial index
	net(:,1) = init(k);
	init(:,k) = []; %remove point k from init to avoid double-counting
	%init = init(:,randperm(length(init))); %optional: random shuffle index

	%%%check each point of init set, add to net if far enough away
	for n = 1:length(init)	
		if far(init(n),net,delta,rho)
			net(:,end+1) = init(n);
		end
	end

	%NOTE: re-write this to be more efficient
	%%%Create struct of neighbors for delta-net: two net points y_n,y_m are connected if
	%rho(y_n,y_m) < 2delta
	N = length(net);
	neighbors = {};
	for n = 1:N
		list =[n]; %first neighbor of y_n will be y_n itself

		%check net for neighbors, add those indices to struct
		%NOTE: this include the point itself as its neighbor
		for m = 1:N
			if rho(net(:,n),net(:,m))<2*delta && m ~= n
				list(end+1) = m;
			end
		end
		neighbors(n).nbr = list;
		
		%TESTING: checking if neighbors are next to net points(they should be)
		%fprintf("neighbors of %f: \n",net(n));
				%for i = 1:length(list)
		%	list(i)
		%	net(:,list(i))
		%end
	end
end

function is_far = far(x,net,delta,rho)
%determine if point 'x' is at least delta-far from set of points 'net'
%inputs: x is vector in R^D, net is D X N matrix, N vectors in net
%output: boolean, true if point is far from set
	is_far = true;
	N = size(net,2);
	k = 1;

	%Check is x is far away from net
	while (is_far && k <= N)
		if rho(x,net(k)) < delta
			is_far = false; %break!
		else
			k = k+1;
		end
	end

end
