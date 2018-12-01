%%%%Main file for implementing the ATLAS algorithm

%S - Simulator, deafult: 1d maggioni example, see simulator.m
%rho - distance function corresponding to output of simulator
%d - dimension of manifold

%Parameters:
%delta - homogenization scale, affects density of sample net
%t_0 - simulation time for short paths (default t_0 = delta^2)
%m - num landmarks for each landmark(m >= d, should be O(d))
%p - num sample paths for each point in net(should be O(delta^-4))
%dt - time step for atlas (default dt = delta/5, should be
%O(delta/ln(1/delta))

d = 1; %for 1-d examples, d = 1
rho = @dist;
delta = 0.1; 
init = [0:0.01:1]; %initial point set for generating delta-net
S = @simulator;
m = 5;
p = 100; %NOTE: p = 10,000 for 1d atlas example, eventually use this
t_0 = 0.01;
%TESTING:
construction(S,init,delta,rho,m,p,t_0,d);

%[net,neighbors] = delta_net(init,delta,rho);
%A = create_landmarks(S,net,m,t_0);
%
%%Plot landmarks for point n
%figure;
%hold on;
%for n  = 1:size(net,2)
%
%	%pPlot each net-point with its landmarks
%	%scatter(A(:,1,n), 0, 'b');
%	%scatter(A(:,2:m+1,n),0*ones(1,m),'r');
%	
%	%optional, 1-d case: plot net-point and landmarks at different heights,
%	%to illustrate that each set of landmarks stays near its net point
%	scatter(A(:,1,n), n, 'b');
%	scatter(A(:,2:m+1,n),n*ones(1,m),'r');
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: start here
function construction(S,init,delta,rho,m,p,t_0,d)
%TODO: comments here
%Currently, this block of code generates delta-net,landmarks, and stores the
%union of the landmarks in a struct array, where each struct index corresponds
%to all neighboring landmarks to a point
	
%Create delta_net from initial points, create landmarks for delta_net
	[net,neighbors] = delta_net(init,delta,rho);
	A = create_landmarks(S,net,m,t_0);

	D = size(net,1); %assuming each column of delta_net is a data_point
	N = size(net,2); %number of points in delta_net
	X = zeros(D,p,N);
	nbr_landmarks = {}; %currently using a struct for collecting neighbor landmarks
	%NOTE: should join this struct with the neighbors struct!
	for n = 1:N
		%simulate p paths around net point y_n
		X(:,1:p,n) = S(net(:,n),p,t_0);

		%NOTE: redundancy here, over collecting landmarks?
		%Take union of all neighboring landmarks to y_n
		net_nbr = neighbors(n).nbr;
		num_nbr = length(net_nbr);

		%store  points from neighboring landmarks as collection of columns
		%associated to y_n
		nbr_landmarks(n).landmark = reshape(A(:,:,net_nbr),D,(m+1)*num_nbr);

		[reduced_landmarks,reduced_X] = LMDS(nbr_landmarks(n).landmark,X,rho,d);

		%Compute new simulator....
		%....

	end
end

function landmarks = create_landmarks(S,delta_net,m,t_0)
%for each y in delta_net, run simulator m times for time t_0 each, keep
%path endpoints as landmarks corresponding to y
%inputs: S - SDE simulator
%		 delta-net - delta-net generated from delta_net function
%		 m - desired number of landmarks to generate for each net point
%		 t_0 - desired time of simulation for each call to S
%outputs: landmarks - array of land-marks, dim D x (m+1) x N
%			1st. dim is coordinate %of landmark point
%			2nd. dim is index of landmark point for given delta-net point x_n,
%		    	index 1 is the delta-net point x_n, next m are landmarks
%			3rd. dim is index of delta_net point
%e.g:
%A(:,:,n) is set of landmarks for delta-net point x_n
%A(:,1,:) is set of delta-net points
	D = size(delta_net,1); %assuming each column of delta_net is a data_point
	N = size(delta_net,2); %number of points in delta_net
	landmarks = zeros(D,m+1,N);
	for n = 1:N
	%generate m paths, save endpoints for delta_net point n
	%NOTE: %first point in landmark list is the %delta_net point itself
			landmarks(:,1,n) = delta_net(n);
			landmarks(:,2:m+1,n) = S(delta_net(n),m,t_0);
		end
	end

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
	%k = 1; %first net point index
	k = randi(length(init)); %optional: random initial index
	net(:,1) = init(k);
	init(:,k) = []; %remove point k from init to avoid double-counting
	init = init(:,randperm(length(init))); %optional: random shuffle index

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
		list =[];	
		%check net for neighbors, add those indices to struct
		for m = 1:N
			if rho(net(:,n),net(:,m))<2*delta
				list(end+1) = m;
			end
		end
		neighbors(n).nbr = list;
	end
end

%
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
function distance = dist(x,y)
%%%Distance function for euclidean space
%inputs: n-dim colummn vectors x,y
	distance = norm(x - y);
end
