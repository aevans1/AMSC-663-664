%%%%TESTING: binning for original simulator
%fprintf("temporarily reducing number of calls to learned simulator \n");
%%%optional parameters
if ~exist('plot_sim','var'), plot_sim = false; end
if ~exist('plot_points','var'), plot_points = false; end
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test: Overall histogram binning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Temporary: setting seed");
seed = 1;
rng(seed);
p = 10000;
T = 50;

%Two ways: random sample initial locations, or all net points
num_locations = 10;

inverse_collection = zeros(d,num_locations*p);
X_collection = zeros(d,num_locations*p);
for m = 1:num_locations
	
	%Uniform sample from [-0.5,1.5];	
	Xzero = 2*rand(d,1) - 0.5;


	%%%Simulate Y paths, bin

	%Find the initial start to start reduced simulator
	for n = 1:size(net,2)
	   	distances(n) = norm(Xzero - net(:,n));
	end
	init_chart = find(distances == min(distances),1);
	
	%Embed initial condition for original simulator so it can berun in reduced
	%simulator
	[local_L,Yzero] = LMDS(L(init_chart).L,Xzero,rho,d);
	
	%Center projected point around the chart center
	Yzero = Yzero - local_L(:,1);
	
	[Ypaths,charts] =learned_simulator(Yzero,p,dt,T,new_S,neighbors,d,delta,net,init_chart);
	
	%Method 1: collect only chart locations
	inverse_Ypaths = [];
	for i = 1:p
		inverse_Ypaths(:,i) = net(:,charts(i));
	end
	inverse_collection(1 + (m-1)*p: m*p) = inverse_Ypaths;
	
	%%%Simulate X paths, bin
	
	init_x = Xzero;
	%simulate p paths around net point y_n
	X = S(init_x,p,T);
	%%bin the endpoints into delta net
	%Method 1: bin the X values by chart
	X_delta = [];
	for i = 1:p
		
		for n = 1:size(net,2)
	    	distances(n) = norm(X(:,i) - net(:,n));
		end
		j =find(distances == min(distances),1);
		X_delta(:,i) = net(:,j);
	end
	X_collection(1 + (m-1)*p: m*p) = X_delta;	
	fprintf("Finished binning samples for init point %d of %d \n",m,num_locations);

end

edges = [-0.5:delta:1.5];

%centers = sort(net);
%d = diff(centers)/2;
%%edges = [centers(1) - d(1) centers(1:end-1) + d centers(end) + d(end)];


figure;
mu_hat = histogram(inverse_collection,edges);
hold on;

mu = histogram(X_collection,edges);
legend('reduced_sim','original_sim');
filename =[datestr(now, 'dd_mmm_yyyy_HH_MM'),'_','d',num2str(d),'_','longtime_test'];
save(filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Test: Individual chart simulators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%N = floor(T/dt);
%
%%path_end will store m column vectors, same dimension as Yzero,
%%collection of endpoints of all simulated paths
%path_end = zeros(d,m); 	
%
%Yzero = rand();
%Y = zeros(d,N+1); 
%Y(:,1) = Yzero;
%
%
%
%%%%Read in chart center for locating points in reduced space
%C = new_S.C;
% for k = 1:m
%        %%%Set initial chart of simulation
%        chart = init_chart;
%		%%%Compute SDE until time dt*N
%		for j = 1:N
%          Y(:,j+1) = learned_simulator_step(Y(:,j),chart,new_S,neighbors,d,dt,delta);
%		end
%		
%		%%%store terminal point in path_end, and store the last chart visited
%		%by that path
%		path_end(:,k) = Y(:,N+1); 
%end
%	if plot_sim	
%			plot_Y( [0:dt:T],Y);
%		end
%	end
%	if plot_points
%		plot_path_end(path_end);
%	end


