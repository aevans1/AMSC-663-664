function gswitch_driver()
%%%%Main file for setting ATLAS parameters

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
	
d = 2;	
rho = @(x,y) norm(x - y);
delta = 7.5; 
m = 10;
p = 10000;
t0 = 0.04;

%dt = t0/5;
%dt = delta/5;
dt = 0.01;

%Set up parameters for original simulator
f = @(x) switch_drift(x);
dt_sim = 0.01; %timestep for original simulator
S = @(Xzero,m,T) simulator(Xzero,m,T,f,dt_sim);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%TESTING: using noisy original system in the beginning construction
% fprintf("%TESTING: using noisy original system in the beginning construction \n");
noise = 1;
% S = @(Xzero,m,T) noisy_simulator(Xzero,m,T,f,dt_sim,noise);
params.noise = noise;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%Data Preprocessing here%%%%%%%%%%%%%%%%%%%
params.net_info = load('gswitch_delta_net.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%save file for use with construction.m

%current_*.mat always contiains most recent _* file run, e.g. driver, atlas
%below lines copy all current outputs into 'gswitch' mat files for use with
%gswitch specific files

save('current_driver.mat','params');
copyfile('current_driver.mat','gswitch_driver.mat');
% construction();
copyfile('current_atlas.mat','gswitch_atlas.mat');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function distance = dist(x,y)
%%%Distance function for euclidean space
%inputs: n-dim colummn vectors x,y
	distance = norm(x - y);
end


function dX = switch_drift(x)
%%%Parameters
k0 = 1; %DNA activation
k1 = 0.0002; %rate of dimerization
gamma_1 = 2;  %rate of de-dimizeration
gamma_0 = 50; %DNA inactivation
gamma_m = 10; %mRNA decay
gamma_n = 1; %protein decay
a = 400;	%transcription
a0 = 0.4;
b = 40;		%translation


%x is a vector in R^3
m = x(1);
n = x(2);
d = x(3);

%
dm = (a0*gamma_0 + a*k0*d)/(gamma_0 + k0*d) - gamma_m*m;
dn = b*m - gamma_n*n - 2*k1*n^2 + 2*gamma_1*d;
dd = k1*n^2 - gamma_1*d;

dX = [dm;dn;dd];

end


