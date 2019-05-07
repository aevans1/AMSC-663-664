function [x,j] = learned_simulator_step(x,i,new_S,neighbors,deg,d,dt,delta)
%Given the learned simulator new_S, computes timestep at point x in chart i
%inputs:    x- initial point for simulator
%           i- chart index that x resides in

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
%           dt = desired time-step length
%           delta- homogenization scale, affects density of sample net
%outputs:   x- new x value after timestep
%           j- new chart index for x

%read in Struct new_S
c = new_S.c;
b = new_S.b;
sigma = new_S.sigma;
T = new_S.T;
mu = new_S.mu;

num_nbr = deg(i);
net_nbr = [i neighbors(i,1:num_nbr)];


%%%Find closest chart to y in chart i coords	
C = squeeze(c(:,i,net_nbr)); %collect centers for chart i
[~,min_dist] = min(sum((x - C).^2,1));
j = net_nbr(min_dist);

if j ~= i
	%%%Set new chart to j
	mu_ij = squeeze(mu(:,i,j));
	mu_ji = squeeze(mu(:,j,i));
	T_ij = squeeze(T(:,:,i,j));
	
	x = T_ij*(x - mu_ij) + mu_ji;
end

%%%Forward Euler step
eta = randn(d,1);
x = x + b(:,j)*dt + sigma(:,:,j)*eta*sqrt(dt);

%%%prevent escape from local chart by applying wall function
if norm(x) > (3/2)*delta
	x = (1/norm(x))*x*(2*delta - (1/2)*delta*exp(3 - (2/delta)*norm(x)));
end

end


