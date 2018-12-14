#ATLAS readme : Work in Progress!

##Main File: atlas.m Description

%%%%Main file for implementing the ATLAS algorithm

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



function [x,j] = learned_simulator_step(x,i,new_S,neighbors,d,dt,delta)
	%TODO: comments here

function [new_Sim,neighbors,net] = construction(S,init,delta,rho,m,p,t_0,d)
%TODO: comments here
	
%Create delta_net from initial points, create landmarks for delta_net
	[net,neighbors] = delta_net(init,delta,rho);

  for n = 1:N
		
    %simulate p paths around net point y_n
		
    
	end

	for n = 1:N
	%construct_SDE(n,neighbors,C,B,Sigma,local_L,local_X,embed_L,embed_X,m,p,t_0)
	%%Compute new simulator....
	
  end
	
		%%%Compute Diffusion coefficients, drift coefficients around y_n
	
    %%%Compute switching maps
		end
	end
	
	%%%Simulator now computed!	
	
	%Now: we have chart centers C, drift coef B, diffusion coeff Sigma
	new_Sim.T = T;
	new_Sim.B = B;
	new_Sim.C = C;
	new_Sim.Sigma = Sigma;
	new_Sim.mu = mu;
	new_Sim.Phi = Phi;	

end

function landmarks = create_landmarks(S,delta_net,m,t_0)
%for each y in delta_net, run simulator m times for time t_0 each, keep
%path endpoints as landmarks corresponding to y
%inputs: S - SDE simulator
%		 delta-net - delta-net generated from delta_net function
%		 m - desired number of landmarks to generate for each net point
%		 t_0 - desired time of simulation for each call to S
%outputs: landmarks - array of land-marks, dim D x (m+1) x N
%			1st. dim is coordinate %of landmark point
%			2nd. dim is index of landmark point for given delta-net point x_n,
%		    	index 1 is the delta-net point x_n, next m are landmarks
%			3rd. dim is index of delta_net point
%e.g:
%A(:,:,n) is set of landmarks for delta-net point x_n
%A(:,1,:) is set of delta-net points

function [net,neighbors] = delta_net(init,delta,rho)
%sub-sample delta-net from given initial set of points, following the
%brute-force method given in section 3.1 of ATLAS paper
%inputs: init - D X N matrix, set of N vectors in R^d%
%		 delta - coarseness of delta-net,
%		 rho -  given distance function
%output: net - columns are data points, all distances >= delta apart
%		 neighbors - struct array, each net point index has a struct with key value
%		 'nbr' that stores the indices of close by net points
%NOTE: all points in domain should be within delta of any net point, not sure
%how to verify?

function is_far = far(x,net,delta,rho)
%determine if point 'x' is at least delta-far from set of points 'net'
%inputs: x is vector in R^D, net is D X N matrix, N vectors in net
%output: boolean, true if point is far from set

function distance = dist(x,y)
%%%Distance function for euclidean space
%inputs: n-dim colummn vectors x,y
	distance = norm(x - y);
end
