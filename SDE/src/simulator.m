function [path_end,regions_visited,X] = simulator(Xzero,m,T,f,dt,regions,dist,plot_sim,plot_points)
%Simulator for linear SDE, using Euler-Maruyama method
%SDE is of form dx = f(x) dt + dW, X(0) = Xzero

%Inputs:
	%Xzero: initial condition for SDE
	%m: number of simulations to run
	%T: right endpoint of time domain
	%f: function taking vector input in R^D, deterministic portion of SDE above
	%dt: timestep 
	%regions: set of points that define regions of state space, 
	%dist: distance to determine if x is close to regions
	%plot_sim: boolean parameter, simulation is plotted
	%plot_path_end: boolean parameter,endpoints of simulations are plotted
%Outputs:	
	%path_end: m column vectors,each the same dimension as Xzero,
		%collection of endpoints of all m simulated paths
	%regions_visited: ith row , j-th col gives region visited by j-th
    	%timestep of i-th path simulated
	
	%%%optional parameters
	if ~exist('plot_sim','var'), plot_sim = false; end
	if ~exist('plot_points','var'), plot_points = false; end
    if ~exist('regions','var'), regions = []; end
    if ~exist('dist','var'), dist = 0 ; regions =[]; end
    %%%

	N = floor(T/dt);

	D = size(Xzero,1);
	path_end = zeros(D,m); 	
	
	X = zeros(D,N+1); 
	X(:,1) = Xzero;

	%%%Calculate m trajectories of SDE
    regions_visited = zeros(1,N+1); %for tracking where paths go in 'regions' set
	for i = 1:m
		%%%Compute SDE until time dt*N
		Xtemp = Xzero;
		for j = 1:N
			dW = sqrt(dt)*randn(D,1); %Brownian random increment in R^D
			Xtemp = Xtemp + f(Xtemp)*dt + dW; %Take Euler-Maryuama step
			X(:,j+1) = Xtemp;
        end
		
        %%%track all of the regions the trajectory visited, if regions were input
        if ~isempty(regions)
           for j = 1:N+1
                regions_visited(i,j) = check_region(X(:,j),regions,dist);
           end
        end        
        
		%%%store terminal point in path_end
        path_end(:,i) = X(:,N+1);
		
		if plot_sim	
			plot_X( [0:dt:T],X);
		end
	end

	if plot_points
		plot_path_end(path_end);
	end

if opts.transition_paths = true
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Original Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute SDE until time dt*N
    %initialize region list, transition tracking
    [dmin,imin] = min(sum((Xzero - regions).^2,1));
    if dmin < dist_sq
        init_region = imin;
    else
        init_region = 0;
    end
    
    original_switches = []; %list of all switches between regions and corresponding times
    switch_count = 0; %tracks number of transitions
    
    tic;
    
    %simulate until certain number of transitions are observed
    while switch_count < max_switches
        t = 0;
        flag = sign(init_region);
        step = 0; %number of steps iterated
        region = init_region;
        
        
        x = Xzero;
     
        %restart simulation once p2 is reached
        while region~=2
            step = step + 1; %number of simulation steps before restart
            
            %Take Euler-Maruyama step
            dW = sqrt(dt_sim)*randn(3,1); %Brownian random increment in R^D
            x = x + switch_drift(x)*dt_sim + sq_noise*dW;
            
            %For now: not tracking charts, just tracking regions
            [~,region_check] = min(sum((x - regions).^2,1));
            
            %only assign a region value if the point is close enough
            if region_check < dist_sq
                new_region = region_check;
            else
                new_region = 0;
            end
            
            if new_region > 0
                if flag == 0
                    region = new_region;
                    t = 0;
                    flag = 1;
                else
                    if region ~= new_region
                        
                        %we've observed a transition!
                        switch_count = switch_count + 1;
                        original_switches = [original_switches;region,new_region,t];
                        t = 0;
                        region = new_region;
                        save('original_tswitch.mat','original_switches','noise');
                        fprintf("transition %d \n",switch_count);
                    end
                end
            end
            t = t+dt_sim;
            
            if mod(step,100000000) == 0
                fprintf("step %d \n",step);
                % 		fprintf("step %d of %f \n",step,N);
                toc;
            end
            
            %increment time-counting step
            
            
        end
        toc;
    end
    original_time = toc;
    %

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%f function for deterministic part, f(x)dt, of SDE above
%function out = f(x)
%    %example, 1d smooth potential
%	%In this case, f(x) = -grad U, where U(x) = 16x^2(x-1)^2
%	out = -32*x*(x-1)*(2*x-1);
%end

%%%plot trajectories of SDE
function plot_X(t,X)
	plot(t,X); hold on;
	xlabel('t');
	ylabel('X');
end

%%%plot terminal points of collection of simulations
function plot_path_end(path_end)
	scatter(zeros(1,length(path_end)),path_end);
end


