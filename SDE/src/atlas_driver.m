function atlas_driver(example)
%%%%Main file for setting ATLAS parameters

%S - Simulator, deafult: 1d maggioni example, see simulator.m
%rho - distance function corresponding to output of simulator
%d - dimension of manifold

%Parameters:
%delta - homogenization scale, affects density of sample net
%t0 - simulation time for short paths (default t0 = delta^2)
%m - num landmarks for each landmark(m >= d, should be O(d))
%p - num sample paths for each point in net(should be O(delta^-4))
%dt - time step for atlas (default dt = delta/5, should be
%O(delta/ln(1/delta))

seed = 1;
fprintf("RNG seed: %d\n",seed);
rng(seed);

if ~exist('example','var') example = 1; end

%Input example determines which system is used for ATLAS
switch example
	
	case 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Example 1: Smooth 1-dimensional Potential from ATLAS paper, ex. 5.2.1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	d = 1;
	rho = @(x,y) norm(x - y);
	delta = 0.1; 
	init = [-0.3:0.01:1.3]; %initial point set for generating delta-net
	
	m = 5;
	p = 10000;
	t0 = delta^2;

	dt = t0/5;
    
	%Set up parameters for original simulator
	f = @(x) example_1_grad(x);
	dt_sim = 0.005; %timestep for original simulator
	S = @(Xzero,m,T) simulator(Xzero,m,T,f,dt_sim);

       
    params.d = d;	
	params.rho = rho;
	params.delta = delta;
	params.m = m;
	params.p = p;
	params.t0 = t0;
	params.dt = dt;
	params.S = S;
	params.dt_sim = dt_sim;
	params.f = f;

    case 2
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Example 2: Rough 1-dimensional Potential from ATLAS paper, ex. 5.2.2
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	d = 1;
	rho = @(x,y) norm(x - y);
	delta = 0.1; 
	init = [-0.3:0.01:1.3]; %initial point set for generating delta-net
	
	m = 5;
	p = 10000;
	t0 = 0.02;
	dt = t0/5;
	
	%Set up parameters for original simulator
	f = @(x) example_2_grad(x);
	dt_sim = 0.00005; %timestep for original simulator
	S = @(Xzero,m,T) simulator(Xzero,m,T,f,dt_sim);

	%As in the paper, run initial point set through simulator for time t= 0.01
	for i = 1:length(init)
		init(:,i) = S(init(:,i),1,0.01);
    end

    params.d = d;	
	params.rho = rho;
	params.delta = delta;
	params.m = m;
	params.p = p;
	params.t0 = t0;
	params.dt = dt;
	params.S = S;
	params.dt_sim = dt_sim;
	params.f = f;

    case 3
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Example 3: Smooth 2-D Potential from ATLAS paper, ex. 5.3.1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	d = 2;	
	rho = @(x,y) norm(x - y);
	delta = 0.2; 
	m = 10;
	p = 10000;
	t0 = delta^2;
	dt = t0/5;

	%Set up parameters for original simulator
	f = @(x) example_3_grad(x);
	dt_sim = 0.005; %timestep for original simulator
    S = @(Xzero,m,T) simulator(Xzero,m,T,f,dt_sim);
	
	%initial  point set for generating delta-net, [-1,2.5] x [-1,2] grid
	%each col of init is (x,y) vector
	init_x = [-1:0.01:2.5];
	init_y = [-1:0.01:2]; 
	[init_X,init_Y] = meshgrid(init_x,init_y);
	temp_init = [init_X(:)';init_Y(:)'];
	Z = [];
	for i = 1:size(temp_init,2)
		Z(i) = example_3(temp_init(:,i));
	end
	
	%Removing all grid points with potential value below 10;
	threshold = [Z < 10; Z < 10];
	init = threshold.*temp_init;
	init = init(:,any(init,1)); %removing zero columns
   
	params.d = d;	
	params.rho = rho;
	params.delta = delta;
	params.m = m;
	params.p = p;
	params.t0 = t0;
	params.dt = dt;
	params.S = S;
	params.dt_sim = dt_sim;
	params.f = f;

    otherwise
	fprintf("enter a number (1-3) \n");
	return;
end

    %Create delta net
    delta_net(init,delta,rho);
	params.net_info = load('current_delta_net.mat');	


fprintf("Running atlas for example %d ... \n ", example);



save('current_driver.mat','params');
construction();



if example == 1 || example == 2
    load('current_atlas.mat');
    load('current_delta_net.mat');
	dim1potential_test(example,new_S,net);
end


%filename =[datestr(now, 'dd_mmm_yyyy_HH_MM'),'_','d',num2str(d),'_','atlas_driver'];
%save(filename);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
%Example 1: 1d smooth potential
%%%%%%%%
function out = example_1(x)
%%
%%example, 1d smooth potential%
    
	out = 16*x^2*(x-1)^2;
end
function out = example_1_grad(x)
    %%% -grad(example_1(x))
	out = -32*x*(x-1)*(2*x-1);
end

%%%%%%%%
%Example 2: 1D rough Potential
%%%%%%%
function out = example_2(x)
	out = 16*x^2*(x-1)^2 + (1/6)*cos(100*pi*x);
end

function out = example_2_grad(x)
    % -grad(example_2(x))
	out = -32*x*(x-1)*(2*x-1) + (1/6)*(100*pi)*sin(100*pi*x);
end

%%%%%%%%
%Example 3: 2d Smooth Potential
%%%%%%%%
function U = example_3(x)
% 2D potential from Crosskey&Maggioni
p1 = [0;0];
p2 = [1.5;0];
p3 = [0.8;1.05];
p = [p1,p2,p3];
c = [0.2,0.2,1/6];
for i = 1 : 3
   d(:,i) = x - p(:,i);
   r(i) = d(:,i)'*d(:,i)/c(i);
   e(i) = exp(-r(i));
end
% potential
U = -log(sum(e));
end

function dU = example_3_grad(x)
% 	%-grad(example_3(x))
% 	%input: x is a vector in R^2
% 	%output: out is a vector in R^2
p1 = [0;0];
p2 = [1.5;0];
p3 = [0.8;1.05];
p = [p1,p2,p3];
c = [0.2,0.2,1/6];
for i = 1 : 3
   d(:,i) = x - p(:,i);
   r(i) = d(:,i)'*d(:,i)/c(i);
   e(i) = exp(-r(i));
end
% potential
% U = -log(sum(e));
s = sum(e);
% gradient
for i = 1 : 3
    grad(:,i) = -2*d(:,i)*e(i)/c(i);
end
dU = (1/s)*(sum(grad,2));
end

function distance = dist(x,y)
%%%Distance function for euclidean space
%inputs: n-dim colummn vectors x,y
	distance = norm(x - y);
end

function dim1potential_test(example,new_S,net) 
%Test function for approximated potentials in 1-dimensional examples 1 and 2 from above.
%Plots approximated potential from ATLAS and compares with true potential
%inputs: example - should be 1 or 2, corresponding to example 1 and example 2
%		new_S - struct that contains ATLAS info, see new_S in construction
%		function
%		net - delta net for ATLAS, see delta_net function
	switch example
		case 1
			U = @(x) 16*x^2*(x-1)^2;
		case 2
			U = @(x) 16*x^2*(x-1)^2 + (1/6)*cos(100*pi*x);
		otherwise 
			fprintf("Please input example = 1 or example = 2 \n");
			return
	end

 	test_set = [-0.25:0.01:1.25];
 	B = new_S.b;
	Sigma = new_S.sigma;
 	Phi = new_S.Phi;
 	
 	%%%effective potential plotting in 1d
 	true_U = [];
 	for k = 1:length(test_set)
 		x = test_set(k);
		true_U(k) = U(x);
 		
		%find closest chart
 		for n = 1:size(net,2)
 			distances(n) = norm(x - net(:,n));
 		end
 		i = find(distances == min(distances),1);
 	
 	    % diffusion coefficent B approximates -grad U, so use -B
		% dividing by Phi, reversing 1d LMDS mapping, i.e. divide by +-1
 		diffs(k) = -B(:,i)/Phi(i);
		sigma(k) = Sigma(:,:,i)/Phi(i);
	end

 	figure

	fprintf("Plotting approximated potential against true potential\n");
	title('potentials');

	n = round(length(test_set)/2);
	piecewise_integrate(0.01,0.5,1,diffs);
 	hold on
	plot(test_set,true_U);
	legend('$\hat{U}(x)$','U(x)');
 	set(legend,'Interpreter','latex');
	
	figure
	title('diffusion coeffs');
	fprintf("Plotting approximated sigma against true sigma \n");
	[sort_net,sort_ind] = sort(net);
	new_sigma = Sigma(:);
	plot(sort_net,new_sigma(sort_ind));
	hold on;
	true_sigma = ones(1,length(test_set));
	plot(test_set,true_sigma);
	ylim([0 1.5]);
	xlim([-1 2]);
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
		if rho(x,net(:,k)) > 0 && rho(x,net(:,k)) < 2*delta 
			is_close = true; %break!
		else
			k = k+1;
		end
	end
end


