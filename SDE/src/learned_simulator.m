function [path_end,chart_end,Y,chart] = learned_simulator(Yzero,m,dt,T,new_S,neighbors,deg,d,delta,net,init_chart,plot_sim,plot_points)
%Euler-Maruyama simulator for 'new_S', the simulator learned from atlas
%inputs:    Yzero- initial point for simulator
%			m- number of simulations to run
%			dt - timestep length
%			T- right endpoint of time domain
%
%           new_S- learned simulator, struct with values:
%				T - d X d X N X N array containing transition mappings between charts,
%				updated each iteration of construct_local_SDE,
%			    T(:,:,i,j) is a d x d matrix which shifts chart coords from chart i to j
%				c - d X N X N array containing local coordinates of neighboring charts,
%				 	updated each iteration of construct_local_SDE,
%				 	c(:,i,j) is a d x 1 vector, isj-th neighbor of net point n,
%		   	  		expressed in chart n coords(LMDS wrt n's landmarks)
%				b - d X 1 vector, drift coordinate for n's local SDE
%				sigma - d x d array, diffusion coords for n's local SDE		
%				mu - d X N X N array containing mean landmarks of local neighboring charts,
%		     		updated each iteration of construct_local_SDE,
%		    		mu(:,i,j) is a d x 1 vector, the average of chart j's landmarks
%		    		expressed in chart n coords
%               Phi- 1-d case only, this is the MDS transformation
%                   	multiplier, Phi(d,j).Phi = +-1, transformation for chart
%                   	y_n

%			neighbors - N x max_deg array, row j is neighbor indices of net	point j
%			deg - N x 1 vector, entry j is number of neighbors(degree) of net point j
%           d- intrinsic dimension of dynamics, problem specifics
%           delta- homogenization scale, affects density of sample net
%			plot_sim- boolean parameter, simulation is plotted
%			plot_points- boolean parameter,endpoints of simulations are plotted
%outputs:   path_end- d x m array of all endpoints from the m simulations
%           chart_end- 1 x m array of all chart locations of path_end
%			Y- d x (T\dt) array of most recent simulations (for m = 1 case)
%			chart- 1 x (T\dt) array of chart locations from Y (for m = 1 case)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%optional parameters

	%plotting
	if ~exist('plot_sim','var'), plot_sim = false; end
	if ~exist('plot_points','var'), plot_points = false; end

	%storing(or not) long trajectories
	trajectory = false; %if true, the trajectory Y will be stored
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	N = floor(T/dt);

	%path_end will store m column vectors, same dimension as Yzero,
	%collection of endpoints of all simulated paths
	path_end = zeros(d,m); 	

	if trajectory
		Y = zeros(d,N+1); 
		Y(:,1) = Yzero;
	else
		%no tracking
		y = Yzero;
		Y = [];
	end

	c = new_S.c;
	b = new_S.b;
	sigma = new_S.sigma;
	T = new_S.T;
	mu = new_S.mu;

	%%%Calculate m trajectories of SDE
    regions_visited = zeros(m,N+1); %for tracking where paths go
	for k = 1:m
        %%%Set initial chart of simulation
        chart(1) = init_chart;
		
		%%%Compute SDE until time dt*N
		for step = 1:N
			
			%%%Take a timestep in the ATLAS learned siulator
			i = chart(step);

			num_nbr = deg(i);
			net_nbr = [i neighbors(i,1:num_nbr)];
			
			
			%%%Find closest chart to y in chart i coords	
			C = squeeze(c(:,i,net_nbr)); %collect centers for chart i
			[~,min_dist] = min(sum((y - C).^2,1));
			j = net_nbr(min_dist);
			
			if j ~= i
				%%%Change coords to new chart j
				mu_ij = squeeze(mu(:,i,j));
				mu_ji = squeeze(mu(:,j,i));
				T_ij = squeeze(T(:,:,i,j));
				
				y = T_ij*(y - mu_ij) + mu_ji;
			end
			
			%%%Forward Euler step
			eta = randn(d,1);
			y = y + b(:,j)*dt + sigma(:,:,j)*eta*sqrt(dt);
			
			%%%prevent escape from local chart by applying wall function
			if norm(y) > (3/2)*delta
				y = (1/norm(y))*y*(2*delta - (1/2)*delta*exp(3 - (2/delta)*norm(y)));
			end
		
		  	chart(step+1) = j;
		
			%only store if path isn't too long
			if trajectory
				Y(:,step+1) = y;
			end

		end
	
		%%%store terminal point in path_end, and store the last chart visited
		%by that path
  		path_end(:,k) = y;
		chart_end(k) = chart(N+1);	
    end
	
	if plot_sim	
			plot_Y( [0:dt:T],Y);
    end

	if plot_points
		plot_path_end(path_end);
	end
	
end

