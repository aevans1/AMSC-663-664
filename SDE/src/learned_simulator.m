function [path_end,chart_end,regions_visited] =learned_simulator(Yzero,m,dt,T,new_S,neighbors,d,delta,net,init_chart,regions,dist,plot_sim,plot_points)
%regions_visited: ith row , j-th col gives region visited by j-th
    %timestep of i-th path simulated
%TODO:comments

%%%optional parameters, for plotting
	if ~exist('plot_sim','var'), plot_sim = false; end
	if ~exist('plot_points','var'), plot_points = false; end
    if ~exist('regions','var'), regions = []; end
    if ~exist('dist','var'), dist = 0 ; regions =[]; end
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
        
        %%%track all of the regions the trajectory visited, if regions were
        %%%input
        if ~isempty(regions)
           for j = 1:N+1
                regions_visited(k,j) = check_region(net(chart(j)),regions,dist);
           end
        end      
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

function region = check_region(x,regions,dist)
   	%inputs: %x: vector in R^D
	%		 %regions: set of values in R^D
	%		 %dist: for checking if x lies near regions, |x - region| < dist?
	%outputs:
		    %region: index of region x lies near, region is 0 if x isn't near
				%any region
    num_regions = size(regions,2);
    found_region = false;
    region = 0; %default is no region
    i = 1;
    while i <= num_regions && ~found_region
        if norm(x - regions(:,i)) < dist
            found_region = true;
            region = i;
        else
            i = i+1;
        end
    end
end





