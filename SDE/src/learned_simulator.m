function [path_end,chart_end] =learned_simulator(Yzero,m,dt,T,new_S,neighbors,d,delta,net,init_chart,plot_sim,plot_points)
%TODO:comments

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

	
    
    %%%Read in chart center for locating points in reduced space
    C = new_S.C;
    
	%%%Calculate m trajectories of SDE
	for k = 1:m
        %%%Set initial chart of simulation
        chart = init_chart;
		%%%Compute SDE until time dt*N
		for j = 1:N
          [Y(:,j+1),new_chart] = learned_simulator_step(Y(:,j),chart,new_S,neighbors,d,dt,delta);
		  chart = new_chart;
		end
		
		%%%store terminal point in path_end, and store the last chart visited
		%by that path
		path_end(:,k) = Y(:,N+1); 
		chart_end(k) = chart;	
	if plot_sim	
			plot_Y( [0:dt:T],Y);
		end
	end

	if plot_points
		plot_path_end(path_end);
	end
end
