function [path_end,chart_end,Y,chart] =learned_simulator(Yzero,m,dt,T,new_S,neighbors,d,delta,net,init_chart,plot_sim,plot_points)
%Euler-Maruyama simulator for 'new_S', the simulator learned from atlas
%inputs:    Yzero- initial point for simulator
%			m- number of simulations to run
%			dt - timestep length
%			T- right endpoint of time domain
%           new_S- learned simulator, struct with values:
%                   T- transition mappings, a struct array
%                   	with T(i,j).T represent the transition matrix from
%                   	chart i to chart j
%                   B- diffusion coefficients for each chart, array of
%                   	form B(:,n) for chart n
%                   C- local coordinates for net points of each chart,
%                   	struct of form C(n,i).C: i-th neighbor of chart point
%                   	y_n, represented in coordinates of chart n
%                   Sigma- drift coefficients for each chart, array of
%                   mu- local coordinates for mean of common landmarks,
%                   	struct array with mu(n,j).mu representing the mean
%                   	landmark for y_j U y_n, represented in chart y_n
%                   Phi- 1-d case only, this is the MDS transformation
%                   	multiplier, Phi(n).Phi = +-1, transformation for chart
%                   	y_n
%           neighbors -struct array, each net point index has a struct with key value
%		             'nbr' that stores the indices of close by net points
%           d- intrinsic dimension of dynamics, problem specifics
%           delta- homogenization scale, affects density of sample net
%i			plot_sim- boolean parameter, simulation is plotted
%			plot_points- boolean parameter,endpoints of simulations are plotted
%outputs:   path_end- d x m array of all endpoints from the m simulations
%           chart_end- 1 x m array of all chart locations of path_end
%			Y- d x (T\dt) array of most recent simulations (for m = 1 case)
%			chart- 1 x (T\dt) array of chart locations from Y (for m = 1 case)

%%%optional parameters, for plotting
	if ~exist('plot_sim','var'), plot_sim = false; end
	if ~exist('plot_points','var'), plot_points = false; end
%%%
	
	N = floor(T/dt);

	%path_end will store m column vectors, same dimension as Yzero,
	%collection of endpoints of all simulated paths
	path_end = zeros(d,m); 	
	
	Y = zeros(d,N+1); 
	Y(:,1) = Yzero;
    
	%%%Calculate m trajectories of SDE
    regions_visited = zeros(m,N+1); %for tracking where paths go
	for k = 1:m
        %%%Set initial chart of simulation
        chart(1) = init_chart;
		%%%Compute SDE until time dt*N
		for j = 1:N
          [Y(:,j+1),new_chart] = learned_simulator_step(Y(:,j),chart(j),new_S,neighbors,d,dt,delta);
		  chart(j+1) = new_chart;
		end
	
		%%%store terminal point in path_end, and store the last chart visited
		%by that path
  		path_end(:,k) = Y(:,N+1); 
		chart_end(k) = chart(N+1);	
    end
	
	if plot_sim	
			plot_Y( [0:dt:T],Y);
    end

	if plot_points
		plot_path_end(path_end);
	end
	
end

%TODO: this may not need to be its own function
function [x,j] = learned_simulator_step(x,i,new_S,neighbors,d,dt,delta)
%Given the learned simulator new_S, computes timestep at point x in chart i
%inputs:    x- initial point for simulator
%           i- chart index that x resides in
%           new_S- learned simulator, struct with values:
%                   T- transition mappings, a struct array
%                   	with T(i,j).T represent the transition matrix from
%                   	chart i to chart j
%                   B- diffusion coefficients for each chart, array of
%                   	form B(:,n) for chart n
%                   C- local coordinates for net points of each chart,
%                   	struct of form C(n,i).C: i-th neighbor of chart point
%                   	y_n, represented in coordinates of chart n
%                   Sigma- drift coefficients for each chart, array of
%                   mu- local coordinates for mean of common landmarks,
%                   	struct array with mu(n,j).mu representing the mean
%                   	landmark for y_j U y_n, represented in chart y_n
%                   Phi- 1-d case only, this is the MDS transformation
%                   	multiplier, Phi(n).Phi = +-1, transformation for chart
%                   	y_n
%           neighbors -struct array, each net point index has a struct with key value
%		             'nbr' that stores the indices of close by net points
%           d - intrinsic dimension of dynamics, problem specifics
%           delta - homogenization scale, affects density of sample net
%           dt = desired time-step length
%outputs:   x- new x value after timestep
%           j- new chart index for x

%read in Struct new_S
C = new_S.C;
B = new_S.B;
Sigma = new_S.Sigma;
T = new_S.T;
mu = new_S.mu;

%TESTING: reading in original chart centers
centers = new_S.centers;

net_nbr = neighbors(i).nbr;
num_nbr = length(net_nbr);

%TESTING: not shifting by centers for no_lmd
sqdist = abs(x - centers(:,net_nbr)).^2;
[~,min_nbr] = min(sum(sqdist,1));
j = net_nbr(min_nbr);
new_chart_center = centers(:,j);

%%%%Code for nearest neighbor with LMDS
%distances = [];
%for n = 1:num_nbr
%    j = net_nbr(n); 
%	distances(n) = norm(x - C(i,j).C);
%end
%[~,min_dist] = min(distances);
%j = net_nbr(min_dist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TESTING: for d = D = 1, no LMDS, with no chart shifting...
%don't need to shift chart coords
%if j ~= i
	%set new chart to i
    %x = (T(i,j).T)*(x - mu(i,j).mu) + mu(j,i).mu;
   
    %TESTING: for d = D = 1 case and no LMDS, changing charts is shifting
    %centers
    %x = x + centers(:,i) - centers(:,j); 
%end

%%%Forward Euler step
eta = randn(d,1);
x = x + B(:,j)*dt + Sigma(:,:,j)*eta*sqrt(dt);

%%%prevent escape from local chart by applying wall function
%TESTING: no_LMDS means need to shift by chart center
%if norm(x) > (3/2)*delta
v = x - new_chart_center;
if norm(v) > (3/2)*delta
	v = (1/norm(v))*v*(2*delta - (1/2)*delta*exp(3 - (2/delta)*norm(v)));
	x = v + new_chart_center;
end
end




