% load('current_delta_net.mat');
% load('current_driver.mat');
load('EM_delta_net.mat');
%%
delta = 5.0;
rho = @(x,y) norm(x - y);
p = 10000;
t0 = 0.5;
dt_sim = 0.01;
f = @(x) drift(x);
S = @(Xzero,m,T) simulator(Xzero,m,T,f,dt_sim);

%%
%%%Read in struct parameters
S = params.S;
delta = params.delta;
rho = params.rho;
m = params.m;
p = params.p;
t0 = params.t0;
d = params.d;


net_info= params.net_info;
net = net_info.net;
neighbors = net_info.neighbors;
edges = net_info.edges;
deg = net_info.deg;
max_deg = net_info.max_deg;
%%
D = size(net,1);
N = size(net,2);
mean_net_distances = zeros(1,N);
t0
p
X = [];
for n = 1:N
    fprintf("Step %d \n of %d", n,N);
    
    %%%Simulate p paths around net point y_n
    X(:,1:p,n) = S(net(:,n),p,t0);
    
    mean_net_distances(n) = diagn_mean_net_dist(n,D,net,delta,t0,X);
end
    
save('mean_net_distances.mat','mean_net_distances');
%%

scatter3(net(1,:),net(2,:),net(3,:),10,mean_net_distances)

%%
function [avg_dist_to_net] = diagn_mean_net_dist(n,D,net,delta,t0,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diagnostic function to check average distances of sample paths from an
%   invidiual net point, and compare with delta value
%Ideally, E[X - net(:,n)] ~ delta
%If this isn't the case, consider changing the sim time t0
%
%inputs:
%		n - current net point, function constructs local_SDE for chart n
%		net - D X N array, columns are data points of delta net
%		t0 - desired time of simulation for each call to S
%       X - D x p array, sample paths from net point(:,n)
%               column i is endpoint of a simulation of time t0
%               starting at delta-net point net(:,n)
%outputs:
%		avg_dist_to_net: approx of E[X - net(:,n)]

%TODO: replace euclidean norm with more general norm

if D == 1
    avg_dist_to_net =  mean(abs(X(:,:,n) - net(:,n)),2);
else
    avg_dist_to_net =  mean(vecnorm(X(:,:,n) - net(:,n)),2);
end
fprintf("avg distance from X to net point n:	%f \n",avg_dist_to_net);
fprintf("	delta: 	%f\n",delta);
fprintf("	t_0: 	%f\n",t0);

end





function dX = drift(x)
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

dm = (a0*gamma_0 + a*k0*d)/(gamma_0 + k0*d) - gamma_m*m;
dn = b*m - gamma_n*n - 2*k1*n^2 + 2*gamma_1*d;
dd = k1*n^2 - gamma_1*d;

dX = [dm;dn;dd];

end
