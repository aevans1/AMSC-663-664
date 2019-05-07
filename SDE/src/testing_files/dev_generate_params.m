
d = 2;	
rho = @(x,y) norm(x - y);
delta = 0.2; 
m = 10;
p = 10000;
t0 = delta^2;
dt = t0/5;

%Set up parameters for original simulator

%f = @(x) example_3_grad(x);
%dt_sim = 0.005; %timestep for original simulator
%S = @(Xzero,m,T) simulator(Xzero,m,T,f,dt_sim);

%initial  point set for generating delta-net, [-1,2.5] x [-1,2] grid
%each col of init is (x,y) vector



%Removing all grid points with potential value below 10;


%delta_net(init,delta,rho);
%params.net_info = load('current_delta_net.mat');	

params.d = d;	
params.rho = rho;
params.delta = delta;
params.m = m;
params.p = p;
params.t0 = t0;
params.dt = dt;
params.S = S;
params.dt_sim = dt_sim;
params.f = f;

return;

