%%%%Main file for implementing the ATLAS algorithm

%S - Simulator 
%rho - distance function corresponding to output of simulator
%d - dimension of manifold

%Parameters:
%delta - homogenization scale, affects density of sample net
%t_0 - simulation time for short paths (default t_0 = delta^2)
%m - num landmarks for each landmark(m >= d, should be O(d))
%p - num sample paths for each point in net(should be O(delta^-4))
%dt - time step for atlass (default dt = delta/5, should be
%O(delta/ln(1/delta))
rho = @dist;
delta = 0.1; %parameter, coarseness of delta net
init = [0:0.01:1]; %initial point set for generating delta-net

function net = delta_net(init,delta,rho)
%%%
%delta_net: sub-sample delta-net from given initial set of points, following the
%brute-force method given in 3.1 of ATLAS papera
%inputs: init - columns are data points,
%		 delta - coarseness of delta-net,
%		 rho -  given distance function
%output: net - columns are data points, all distances >= delta apart
%NOTE: all points in domain should be within delta of any net point, not sure
%how to verify?
%%%
	k = 1; %first net point index
	%k = randi(length(init)); %optional: random initial index
	net(:,1) = init(k);
	init(:,k) = []; %remove point k from init
	init = init(:,randperm(length(init))); %optional: random shuffle index

	%%%check each point of init set, add to net if far enough away
	for n = 1:length(init)	
		if far(init(n),net,delta,rho)
			net(:,end+1) = init(n)
		end
	end
end

%NOTE: start here
function landmarks = create_landmarks(delta_net,m,t_0)
	%for each y in delta_net, run simulator m times for time t_0 each, keep
	%path endpoints as landmarks corresponding to y
	
		landmarks = [];
end



%%%determine if point is at least dellta-far from set
%inputs: x is vector in R^n, net is n X N matrix, N vectors in net
%output: boolean, true if point is far from set
function is_far = far(x,net,delta,rho)
	is_far = true;
	N = size(net,2);
	k = 1;

	%Check is x is far away from net
	while (is_far && k <= N)
		if dist(x,net(k)) < delta
			is_far = false; %break!
		else
			k = k+1;
		end
	end

end
%%%Distance function for euclidean space
%inputs: n-dim colummn vectors x,y
function distance = dist(x,y)
	distance = norm(x - y);
end
