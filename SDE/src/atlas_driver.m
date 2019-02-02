%%%%Main file for setting ATLAS parameters

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

close all;

%seed rng
rng(sum(clock));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 1: Smooth 1-dimensional Potential from ATLAS paper, ex. 5.2.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%d = 1;
%rho = @dist;
%delta = 0.1; 
%dt = delta/5;
%init = [-0.3:0.01:1.3]; %initial point set for generating delta-net
%
%m = 5;
%p = 10000;
%t_0 = 0.01;
%
%%Set up parameters for original simulator
%f = @(x) example_1_grad(x);
%dt_original = 0.005; %timestep for original simulator
%S = @(Xzero,m,T) simulator(Xzero,m,T,f,dt_original);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 2: Rough 1-dimensional Potential from ATLAS paper, ex. 5.2.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%d = 1;
%rho = @dist;
%delta = 0.1; 
%dt = delta/5;
%init = [-0.3:0.01:1.3]; %initial point set for generating delta-net
%S = @simulator;
%m = 5;
%p = 10000;
%t_0 = 0.02;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 3: Smooth 2-D Potential from ATLAS paper, ex. 5.3.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 2;
rho = @dist;
delta = 0.2; 

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
m = 5;
p = 10000;
t_0 = delta^2;
dt = t_0/5;

%Set up parameters for original simulator
f = @(x) example_3_grad(x);
dt_original = 0.005; %timestep for original simulator
S = @(Xzero,m,T) simulator(Xzero,m,T,f,dt_original);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[new_S,neighbors,net] = construction(S,init,delta,rho,m,p,t_0,d);

%Yzero = rand();
%T = 10;
%Ypaths =learned_simulator(Yzero,p,dt,T,new_S,neighbors,d,delta,net);

%figure;
%histogram(Ypaths,10);
%hold on;
%%%%TESTING: binning for original simulator
%% init_x = rand();
%init_x = Yzero;
%%simulate p paths around net point y_n
%X = S(init_x,p,T);
%histogram(X,10);

%%%TESTING 
%%%Construction of Charts and Effective Potential, 1D 
% 
%  [new_S,neighbors,net] = construction(S,init,delta,rho,m,p,t_0,d);
% % 
 %Constructing effective potential
% test_set = [-0.25:0.01:1.25];
% B = new_S.B;
% Phi = new_S.Phi;
% 
% %%%effective potential plotting
% U = [];
% for k = 1:length(test_set)
% 	x = test_set(k);
% 	%find closest chart
% 	for n = 1:size(net,2)
% 		distances(n) = norm(x - net(:,n));
% 	end
% 	i = find(distances == min(distances),1);
% 	%%invert MDS if d = 1
% 	%d(k) = B(i)/Phi(i).Phi;
% 	%B approximates -grad U, so use -B
% 	diffs(k) = -B(:,i)/Phi(i).Phi;
% 	%fprintf("experiment: \n");
% 	%diffs(k) = -t_0*B(:,i)/Phi(i).Phi;
% 	U(k) = 16*x^2*(x-1)^2;
% 	%V(k) = 16*test_set(k)^2*(x-1)^2 + (1/6)*cos(100*pi*test_set(k));
% end
% figure
% 
%%OLD 
% %piecewise_integrate(0.01,-0.25,1.25,diffs);
% 
% %Uncomment for example 2
% piecewise_integrate(0.01,0.5,1,diffs);
% hold on;
% 
% plot(test_set,U);
% 
% %Uncomment for test 2, rough potential
% %plot(test_set,V);
% 
% legend('$\hat{U}(x)$','U(x)');

%Uncomment for test 2, rough potential
%legend('$\hat{V}(x)$','V(x)');

% set(legend,'Interpreter','latex');

%%%%TESTING
%%%%simulator
%x = 0.2;
%Y(1) = x;
% for t = 1:9
	%find closest chart
% 	for n = 1:size(net,2)
% 		distances(n) = norm(Y(t) - net(n));
% 	end
% 	i = find(distances == min(distances));
% 	Y(t+1) = learned_simulator_step(Y(t),i,new_S,neighbors,d,dt,delta);
% end

%plot([0:9],Y);
%%%%

%%%TESTING
%%%delta-net, landmarks
%[net,neighbors] = delta_net(init,delta,rho);
%A = create_landmarks(S,net,m,t_0);
%
%figure(1); hold on;
%
%for n = 1:size(net,2) 
%    scatter(net(:,n), 0, 30,'b','filled');
%end
%
%ylim([-1,size(net,2) + 1]);
%xlim([-0.2,1.2]);
%title(['delta-net on [0,1], delta = 0.1']);
%%Plot landmarks for point n
%figure(2);
%hold on;
%for n  = 1:size(net,2)
%
%	%pPlot each net-point with its landmarks
%	scatter(A(:,1,n), 0,30, 'b','filled');
%	scatter(A(:,2:m+1,n),0*ones(1,m),20,'r','filled');
%	
%	%optional, 1-d case: plot net-point and landmarks at different heights,
%	%to illustrate that each set of landmarks stays near its net point
%	scatter(A(:,1,n), n, 30 ,'b','filled');
%	scatter(A(:,2:m+1,n),n*ones(1,m),20,'r','filled');
%end
%ylim([-1,size(net,2) + 1]);
%xlim([-0.2,1.2]);
%title(['delta-net with landmarks']);
%%%

%%%%%%%%
%Example 1: 1d smooth potential
%%%%%%%%
function out = example_1(x)
    %example, 1d smooth potential
	out = 16*x^2*(x-1)^2;
end
function out = example_1_grad(x)
    % -grad(example_1(x))
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
function out = example_3(x)
	%input: x is a vector in R^2
	%output: out is a scalar
	c = [1/5;1/5;1/6];
	p_1 = [0;0];
	p_2 = [1.5;0];
	p_3 = [0.8;1.05];
	out =  -log( exp( (1/c(1))*(-norm(x - p_1)^2)) + ...
		exp( (1/c(2))*(-norm(x - p_2)^2)) + ...
		exp( (1/c(3))*(-norm(x - p_3)^2)));
end

function out = example_3_grad(x)
	%-grad(example_3(x))
	%input: x is a vector in R^2
	%output: out is a vector in R^2
	
	c = [1/5;1/5;1/6];
	p_1 = [0;0];
	p_2 = [1.5;0];
	p_3 = [0.8;1.05];
	out =  -2*( ( (x - p_1)/c(1))*exp( (1/c(1))*(-norm(x - p_1)^2)) + ...
		( (x - p_2)/c(2))*exp( (1/c(2))*(-norm(x - p_2)^2)) +  ...
		((x-p_3)/c(3))*exp( (1/c(3))*(-norm(x - p_3)^2)))/...
		(exp( (1/c(1))*(-norm(x - p_1)^2)) + ...
		exp( (1/c(2))*(-norm(x - p_2)^2)) + ...
		exp( (1/c(3))*(-norm(x - p_3)^2)));
end

function distance = dist(x,y)
%%%Distance function for euclidean space
%inputs: n-dim colummn vectors x,y
	distance = norm(x - y);
end
