function atlas_driver(example)
rng(sum(clock));
	d = 1;
	rho = @(x,y) norm(x - y);
	delta = 0.15; 
	init = [-10:0.1:10]; %initial point set for generating delta-net
	
	[net_neighbors] = delta_net(init,delta,rho);
