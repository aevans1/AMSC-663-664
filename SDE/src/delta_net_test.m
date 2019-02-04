%shared variable section

%Uncomment one of 2d or 3d example
%2D example
%x_init = [-1:0.01:2.5];
%y_init = [-1:0.01:2];
%[X,Y] = meshgrid(x_init,y_init);
%init = [X(:)' ; Y(:)'];
fprintf("Running test for a 2D grid example \n");

%3D grid example
x_init = [-1:0.1:1];
y_init = x_init;
z_init = x_init;
[X,Y,Z] = meshgrid(x_init,y_init,z_init);
init = [X(:)'; Y(:)'; Z(:)' ];
fprintf("Running test for a 3D grid example \n");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = @(p_1,p_2) norm(p_1 - p_2);
delta = 0.5;
[net,neighbors] = delta_net(init,delta,rho);
N = size(net,2);

plot_net(net,neighbors);
%TODO: write a plot function for delta net and neighbors
%2D Example, plot for fun


%% Test 1: Delta Spacing
%TODO: write fail, pass meanings
copy_net = net;
n = 1;

while (n < N)
	copy_net(:,1)= [];
	delta_spaced = far(net(:,n),copy_net,delta,rho); %is point far away?
	assert(delta_spaced);
	n = n + 1;
end

%% Test 2: Neighbors
%TODO: write fail, pass meanings
n = 1;
while (n < N)
	net_nbr = neighbors(n).nbr; %global indices of neighbors of net point idx n
	organized_neighbors = ( all(net(:,net_nbr(1)) == net(:,n)) );
	assert(organized_neighbors); %1st index nbr should be the point itself
	num_nbr = length(net_nbr);
	for i = 2:num_nbr
		m = net_nbr(i); %global index of net point nbr
		spaced_neighbors = (rho(net(:,n),net(:,m)) < 2*delta);
		assert(spaced_neighbors); %neighbors should be within 2 delta
	end
	n = n + 1;
end

%NOTE; this test doesn't cover the case for where there are NO net points in a
%region of the space
%% Test 3: Delta Covering 
%(depends on neighbor tests, test 2)
%Fail: there is a delta_net point with no non-trivial(not itself) neighbors assigned
%This means: 1) the neighbor construction function is wrong
%			 2) there is no other net point within 2delta of above point,
%			 	and the domain is not covered adequately by the delta_covering,
%			 	user needs to make a new delta_net
%True: every delta_net point has non-trivial neighbors, according to the
%	   neighborhood construction function
delta_cover = true;
n = 1;
while (n < N)
	net_nbr = neighbors(n).nbr;
	num_nbr = length(net_nbr);
	assert(num_nbr > 1);
	n = n + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function is_far = far(x,net,delta,rho)
%%NOTE: cross-ref with far function in delta_net.m
%%NOTE: if input x is part of net, this will return false
	
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

function plot_net(net,neighbors)

	%check for all inputs
	assert(exist('neighbors','var'),'This function needs a neighbor struct');
		
	D = size(net,1); %dimension of vectors in the net
	switch D
			case 1
				plotter = 'plot( [plot_array(1,:)],''-o'',''Color'',''k'');';
			case 2
				plotter = ['plot( [plot_array(1,:)],[plot_array(2,:)],'...
				'''-o'',''Color'',''k'');'];
			case 3
				plotter = ['plot3( [plot_array(1,:)],[plot_array(2,:)],'...
				'[plot_array(3,:)],''-o'',''Color'',''k'');'];
	end	
	for n = 1:size(net,2)
		nbrs = neighbors(n).nbr;
		num_nbr = size(nbrs,2);
		for i = 1:num_nbr
			nbr = nbrs(i);
			plot_array = net(:,[n nbr]);
			eval(plotter);
			hold on;
		end
	end
end
