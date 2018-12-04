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
dt = delta/5;
init = [0:0.01:1]; %initial point set for generating delta-net
S = @simulator;
m = 5;
p = 100; %NOTE: p = 10,000 for 1d atlas example, eventually use this
t_0 = 0.01;
%TESTING:
[new_S,neighbors,net] = construction(S,init,delta,rho,m,p,t_0,d)

test_set = [0:0.03:1]
for k = 1:length(test_set)
	x = test_set(k);
	
	%find closest chart
	for n = 1:size(net,2)
		distances(n) = norm(x - net(n));
	end
	i = find(distances == min(distances));
	
	B = new_S.B;
	U(k) = B(:,i);
end

plot(test_set,U);



%TESTING simulator
%x = 0.2;
%%find closest chart
%for n = 1:size(net,2)
%	distances(n) = norm(x - net(n));
%end
%i = find(distances == min(distances));
%
%Y(1) = x;
%for t = 1:9
%	%find closest chart
%	for n = 1:size(net,2)
%		distances(n) = norm(Y(t) - net(n));
%	end
%	i = find(distances == min(distances));
%	Y(t+1) = learned_simulator_step(Y(t),i,new_S,neighbors,d,dt,delta);
%end
%
%plot([0:9],Y);

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
%NOTE: for now, passig in struct that give neighbor inndices, 'neighbors'
function [x,j] = learned_simulator_step(x,i,new_S,neighbors,d,dt,delta)
	%read in Struct new_S
	C = new_S.C;
	B = new_S.B;
	Sigma = new_S.Sigma;
	T = new_S.T;
	mu = new_S.mu;

	net_nbr = neighbors(i).nbr;
	num_nbr = length(net_nbr);
	distances = [];
	for n = 1:num_nbr
		j = net_nbr(n);
		distances(n) = norm(x - C(n,j).C);
	end	

	%find global index for closest chart
	temp = find(distances == min(distances));
	new_i = net_nbr(temp);

	mu(new_i,i).mu;
	if new_i ~= i
		x = (T(i,new_i).T).' * (x - mu(i,new_i).mu) + mu(new_i,i).mu;
	end

	%Forward Euler step
	%take eta from gaussian N(0,eye(d))
	eta = (randn(d)).';
	x = x + B(:,i)*dt + (eta.')*Sigma(:,i)*sqrt(dt);

	%prevent escape from local chart
	%NOTE: I don't know how to write this for vectors x, just scalars(?)
	if norm(x) > (3/2)*delta
		x = (1/norm(x))*x*(2*delta - (1/2)*delta*exp(3 - (2/delta)*norm(x)));
	end

end


function [new_Sim,neighbors,net] = construction(S,init,delta,rho,m,p,t_0,d)
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
	%nbr_landmarks = {}; %currently using a struct for collecting neighbor landmarks
	%NOTE: should join this struct with the neighbors struct!
	
	
	%TEMP: better way to write this?
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
		%nbr_landmarks(n).L = reshape(A(:,:,net_nbr),D,(m+1)*num_nbr);
		L = reshape(A(:,:,net_nbr),D,(m+1)*num_nbr);

		[embed_L(n).L,embed_X(n).X] = LMDS(L,X,rho,d);
		fprintf("original: \n");
		L
		fprintf("new: \n");
		embed_L(n).L	
	end

	for n = 1:N
	%Compute new simulator....
		net_nbr = neighbors(n).nbr;
		num_nbr = length(net_nbr);
		
		local_L = embed_L(n).L;
		local_X = embed_X(n).X;

		%select neighboring delta-net points:
		%local_L = [y_n |y_n landmarks |nbr1 |nbr 1landmarks |nbr2|nbr2
		%landmarks...];
		ind = [1:m+1:(m+1)*num_nbr];

		%%%Compute chart 'C_n' for y_n and center around embedded y_n
		testing = [];
		for i = 1:num_nbr
			j = net_nbr(i); %global index of neighbor
			C(n,j).C = local_L(:,ind(i)) - local_L(:,1); %find nbr and center
		end
		
		%%%Compute Diffusion coefficients, drift coefficients around y_n
		B(:,n) = (1/p*t_0)*sum(local_X,2); %drift for y_n
		Sigma(:,n) = sqrt((1/t_0)*cov(local_X.',1)); %diffusion for y_n
		
		%%%Compute switching maps
		%TODO: re-write this linea
		%Grab max index of nbr with global index less than n
		%Switch ind corresponds to indices of L matrix rather than delta_net
		
		%switch_ind = (net_nbr < n).*[1:num_nbr];
		%switch_ind = max(switch_ind);
		switch_ind = net_nbr(net_nbr < n);
	
		%NOTE: net_nbr = [n idx1 idx 2 ...], idx1 < idx 2 < ...

		%NOTE: change this, too confusing, but switch_ind+1 is index of y_n
		%landmarks in L matrix

		%grab all landmarks for neighbors with index less than L
		%nbr_L = local_L(:,1:switch_ind*(m+1)); 
		%NOTE: change in,ni notation, this is terrible
		for i = 1:length(switch_ind)
			
			j = switch_ind(i);

			local_nbr_L = embed_L(j).L; %embeddings wrt nbr of y_n

			%NOTE: indexing to find nbr is based on that net_nbr = [n idx1 idx2
			%.] and switch_ind = [idx1 idx1 ...], switch_ind + 1 will be index
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
			%NOTE:*pseudo inv,  (X^*Y)^* = Y^*X
			%T(j,n).T = pinv(L_jn - mu_jn)*(L_nj - mu_nj) = pinv(T(n,j).T)


			%TODO: check linear algebra here
			T(n,j).T = pinv( (L_nj - mu(n,j).mu).' )*( (L_jn - mu(j,n).mu).' );
			T(j,n).T = pinv(T(n,j).T);
		end
	end
	%%%Simulator now computed!	
	%TESTING
	%for n = 1:N
	%	fprintf("neighbors of %d: \n",n);
	%	nbr = neighbors(n).nbr
	%	fprintf("transitions from index %d \n",n);
	%	for i = 1:length(nbr)
	%		fprintf("%d to %d \n",n,nbr(i));
	%		T(n,nbr(i)).T	
	%	end
	%end
	
	%Now: have chart centers C, drift coef B, diffusion coeff Sigma
	new_Sim.T = T;
	new_Sim.B = B;
	new_Sim.C = C;
	new_Sim.Sigma = Sigma;
	new_Sim.mu = mu;
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
		list =[n]; %first neighbor of y_n will be y_n itself

		%check net for neighbors, add those indices to struct
		%NOTE: this include the point itself as its neighbor
		for m = 1:N
			if rho(net(:,n),net(:,m))<2*delta && m ~= n
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
