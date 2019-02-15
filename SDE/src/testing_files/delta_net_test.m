function delta_net_test(net_params)
%Contains tests for an input delta net to assure that it meets the
%qualifications for a delta net, and that the neighbors of the delta net are
%sufficiently close
%inputs: a struct net_params with key values
%		net
%		neighbors
%		delta
%		rho	
%	where all the above are the same parameters from a construction of the
%	function [net,neighbors] = delta_net(init,delta,rho)
%	If no user input, default is to use a 2D delta net example and test that
%	net
%outputs: assert statements regarding three tests described below
%TODO: add test for condition that delta_net covers the initial space	
	if ~exist('net_params','var') 
		%%%%If no user input, generate a net from an example	
		x_init = [-1:0.01:2.5];
		y_init = [-1:0.01:2];
		[X,Y] = meshgrid(x_init,y_init);
		init = [X(:)' ; Y(:)'];
		fprintf("Running test for a 2D grid example \n");

		rho = @(p_1,p_2) norm(p_1 - p_2);
		delta = 0.5;
		[net,neighbors] = delta_net(init,delta,rho);
		
		%save the above for later, just in case
		net_params.net = net;
		net_params.neighbors = neighbors;
		net_params.rho = rho;
		net_params.delta = delta;
		save('delta_net_test');
	else
		net = net_params.net;
		neighbors = net_params.neighbors;
		rho = net_params.rho;
		delta = net_params.delta;
	fprintf("Running test for user input delta_net \n");

	end	
	
	N = size(net,2);
	
	plot_net(net,neighbors);
	%TODO: write a plot function for delta net and neighbors
	%2D Example, plot for fun
	
	
	%% Test 1: Delta Spacing
	%Pass: all of the net points are at least delta away from each other
	%Fail: there is at least one pair of points closer than delta
	copy_net = net;
	n = 1;
	
	while (n < N)
		copy_net(:,1)= [];
		delta_spaced = far(net(:,n),copy_net,delta,rho); %is point far away?
		try 
			assert(delta_spaced);
		catch
			fprintf(['Error 1: net point %d is too close to another net point \n'],n);
		end
		%assert(delta_spaced,msg);
		n = n + 1;
	end
	
	%% Test 2: Neighbors
	%Pass: Passed both asserts for all net points
	%		assert 'organized neighbors' pass for point n:
	%			the array of neighbors of net point n has first column as net
	%			point n (this structure is assumed in later functions and is
	%			needed)
	%		assert 'spaced neighbors' pass for point n:
	%			all neighbors in the array of neighbors for point n are within
	%			2*delta
	%Fail: There is a least one net point such that doesn't pass both asserts
	%above, details in error statement
	n = 1;
	while (n < N)
		net_nbr = neighbors(n).nbr; %global indices of neighbors of net point idx n
		organized_neighbors = ( all(net(:,net_nbr(1)) == net(:,n)) );
		try
			assert(organized_neighbors);
		catch
			fprintf(['Error 2: net point %d is not listed as its first '...
			'neighbor \n'],n);
		end
		num_nbr = length(net_nbr);
		for i = 2:num_nbr
			m = net_nbr(i); %global index of net point nbr
			spaced_neighbors = (rho(net(:,n),net(:,m)) < 2*delta);
			try
				assert(spaced_neighbors);
			catch
				fprintf(['Error 3: net points %d and %d are neighbors but are too'...
				' far \n'],n,m);
			end
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
			otherwise
				fprintf("Dimension Error: this function only plots 1D, 2D, or 3D data \n");
				return
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
