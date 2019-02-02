%Simulator for linear SDE, using Euler-Maruyama method
%SDE is of form dx = f(x) dt + dW, X(0) = Xzero
	
	%mu: controls random portion of sde above
	%f: function taking one input of real numbers, deterministic portion of SDE above, e.g f(x) = -32*x(x-1)*(2x-1)
	%Xzero: initial condition for SDE, in [0,1], e.g. 0.3
	%T: right endpoint of time domain [0,T], e.g. 50
	%N: number of timesteps, e.g: 2^15
	%det: boolean, 1 means compute deterministic DE, 0 means don't
	%Xem: Euler-Maryuma solution of SDE above
	%X: Eulers method solution of DE dx/dt = f(x)
	%plot_sim: boolean parameter, simulation is plotted
	%plot_path_end: boolean parameter,endpoints of simulations are plotted
	
%TODO: comments for simulator 
function [path_end] = simulator(Xzero,m,T,f,dt,plot_sim,plot_points)

	%%%optional parameters
	if ~exist('plot_sim','var'), plot_sim = false; end
	if ~exist('plot_points','var'), plot_points = false; end
	%%%

	N = floor(T/dt);

	D = size(Xzero,1);
	%path_end will store m column vectors, same dimension as Xzero,
	%collection of endpoints of all simulated paths
	path_end = zeros(D,m); 	
	
	X = zeros(D,N+1); 
	X(:,1) = Xzero;

	%%%Calculate m trajectories of SDE
	for i = 1:m

		%%%Compute SDE until time dt*N
		Xtemp = Xzero;
		for j = 1:N
			dW = sqrt(dt)*randn(D,1); %Brownian random increment in R^D
			Xtemp = Xtemp + f(Xtemp)*dt + dW; %Take Euler-Maryuama step
			X(:,j+1) = Xtemp;
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



