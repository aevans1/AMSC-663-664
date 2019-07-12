function [net_info] = delta_net(init,delta,rho,is_random)

%sub-samples initial set of points to create a delta-net, following the
%brute-force method given in section 3.1 of ATLAS paper
%inputs: init - D X N matrix, set of N vectors in R^d%
%		 delta - coarseness of delta-net,
%		 rho -  given distance function
%output: net_info - struct for delta_net, contains
%			(N = number of net points)
%			net - d X N array, columns are data points of delta net
%			neighbors - N x max_deg array, row j is neighbor indices of net	point j
%			deg - N x 1 vector, entry j is number of neighbors(degree) of net point j
%			max_deg - maximum enty of deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optional parameters

%%%Default mode is random
	if ~exist('is_random','var'), is_random = true; end

	%%%Random mode shuffles initial point set
	if is_random
		init = init(:,randperm(length(init)));
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%Step 1: Build the Delta Net
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
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

	%Step 2: Find and store neighbors of each delta net point
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%Create network from delta-net:
	%%%two net points y_i,y_j are connected if rho(y_i,y_j) < 2delta
	figure; grid; hold on;
	N = size(net,2);
    D = size(net,1);
	
    edges = []; %edges of the delta-net
	for i = 1 : N - 1
		for j = i + 1 : N
			if rho(net(:,i),net(:,j)) < 2*delta
				edges = [edges;[i j]];
               
                if D == 2
                    plot(net(1,[i,j]),net(2,[i,j]),'color','b','Linewidth',0.5);
                elseif D == 3
                    plot3(net(1,[i,j]),net(2,[i,j]),net(3,[i,j]),'color','b','Linewidth',0.5);
                end
                
			end
		end
	end

	%%%Find degrees of each net point, maximal degree of net, remove isolated points
	deg = zeros(N,1);
	num_iso = 0;
	isolated = [];
	for n = 1 : N
		ind = find(edges(:,1) == n | edges(:,2) == n);
		if isempty(ind)
			num_iso = num_iso + 1;	
			fprintf("removing isolated point!consider re-initializing delta_net \n");
			isolated(num_iso) = n;
		else
			deg(n) = length(ind);
		end
	end
	max_deg = max(deg);

	if ~isempty(isolated)
		net(:,isolated) = [];
		deg(isolated) = [];
		N = N - num_iso;
	end

	%%%Create neighbor set from edge set
	neighbors = zeros(N,max_deg);
	for n = 1 : N
		ind = find(edges(:,1) == n);
		if ~isempty(ind)
			l1 = length(ind);
			neighbors(n,1 : l1) = edges(ind,2)';
		else
			l1 = 0;
		end

		ind = find(edges(:,2) == n);
		if ~isempty(ind)
			l2 = length(ind);
			neighbors(n, l1 + 1 : l1 + l2) = edges(ind,1).';
		end
    end

    %set good view for 3D
    if (D == 3)
        view(3); 
    end
    
 	net_info.net = net;
 	net_info.neighbors = neighbors;
 	net_info.edges = edges;
    net_info.deg = deg;
	net_info.max_deg = max_deg;
    

	save('current_delta_net.mat','net','neighbors','edges','deg','max_deg');
end

function is_far = far(x,net,delta,rho)
%%%NOTE: if input x is part of net, this will return false

%determines if point 'x' is at least delta-far from set of points 'net'
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

