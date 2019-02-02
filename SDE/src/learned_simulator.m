function [path_end] =learned_simulator(Yzero,m,dt,T,new_S,neighbors,d,delta,net,plot_sim,plot_points)
%TODO:comments

%%%optional parameters, for plotting
	if ~exist('plot_sim','var'), plot_sim = false; end
	if ~exist('plot_points','var'), plot_points = false; end
	%%%
	
	N = floor(T/dt);

	%path_end will store m column vectors, same dimension as Yzero,
	%collection of endpoints of all simulated paths
	path_end = zeros(size(Yzero,1),m); 	
	
	Y = zeros(1,N+1); 
	Y(1) = Yzero;

    
    %%%Read in chart center for locating points in reduced space
    C = new_S.C;

    
	%%%Calculate m trajectories of SDE
	for k = 1:m

        %%%Set initial chart of simulation
        %%%TODO: figure out how to do this rigorously?
        i = 1;
		%%%Compute SDE until time dt*N
		for j = 1:N
          [Y(j+1),i] = learned_simulator_step(Y(j),i,new_S,neighbors,d,dt,delta);
		end
		
		%%%store terminal point in path_end
		path_end(:,k) = Y(N+1); 
		
	if plot_sim	
			plot_Y( [0:dt:T],Y);
		end
	end

	if plot_points
		plot_path_end(path_end);
	end
end
