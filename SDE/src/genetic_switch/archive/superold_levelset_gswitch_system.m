%%%Parameters and Equations take from Lv et al 2014 Paper
%%%Simulation in 3D of perturbed ODE for the mean field ODE satisfied by mRNA,
%proteins and dimers
dt = 0.01;
D = 3; %dimension of problem
epsilon = 20; %noise parameter
num_points = 1; %number of trajectories for each sim. starting point

%%%%These lines are for if a set size net is desired
%num_net = 50; %number of delta net points desired
%net = zeros(D,num_net);
net = [];

%rho = @(x,y) ellipse_dist3(x,y,1.0,80.0,4.0); %distance function for delta-net
rho = @(x,y) norm(x - y);
delta = 0.5; %spatial homogenization parameter, delta-net

scale = [1; 80; 4];

%xi - asympt. stable equilibria representing inactive switch
%xa - asympt. stable equilibria representing active switch
%xs - morse index one saddle
xi = [0.0402067142317042; 1.60826856926817; 0.000258652779089588];
xa = [29.3768600805981; 1175.07440322392; 138.079985311206];
xs = [10.5829; 423.3173; 17.9198];

scaled_xi = xi./scale;
scaled_xa = xa./scale;
scaled_xs = xs./scale;

figure;
hold on;
view(3); grid;

%%%Plot re-scaled saddle and equilibria
plot3(xi(1),xi(2),xi(3),'.','Color','r','MarkerSize',40);
plot3(xa(1),xa(2),xa(3),'.','Color','g','MarkerSize',40);
plot3(xs(1),xs(2),xs(3),'.','Color','y','MarkerSize',40);

%load in initial point sets for generating delta net
data = load('reduced_QPotLevelSets.mat');
reduced_patches = data.reduced_patches;
col = [winter(3);summer(3)];

for k = 1 : 6
    
    %read in path, plot
    my_field = strcat('v0_',num2str(k));
    v = reduced_patches.(my_field);
    my_field = strcat('f0_',num2str(k));
    f = reduced_patches.(my_field);
    %p = patch('Vertices',v,'Faces',f,'Facecolor',col(k,:),'Edgecolor',col(k,:));
    alpha(0.4);
    
    %scale patch, start simulating
    scaled_v = (v')./scale;
    num_v = size(scaled_v,2);
    
    %%%Adjust num of trajectories depending on the level set
    %require more exploration for patches 4 - 6, they're sparser patches
    if k <= 3
        num_points = 1;
    else
        %num_points = 5;
        num_points = 1;
    end
    
    for n = 1 : num_v
        x = scaled_v(:,n);
        
        %%%Compute SDE until time dt*num_steps
        N = 0;
        while N < num_points
            %Euler-Maruyama step
            dW = sqrt(dt)*randn(D,1); %Brownian random increment in R^D
            x = x + scaled_switchgrad(x)*dt + sqrt(epsilon)*dW;
            
            %Add point to net if delta away from net
            if far(x,net,delta,rho)
                N = N + 1; %update length of net
                net = [net x];
            end
        end
    end
    fprintf("level set %d finished \n",k);
end
%%

rescaled_net = net.*scale;
plot3(rescaled_net(1,:),rescaled_net(2,:),rescaled_net(3,:),'.');

params.dt = dt;
params.epsilon = epsilon;
params.delta = delta;
params.rho = rho;
params.num_net = N; %stored num_net is furthest step the loop got to

%save('gswitch_data','net','rescaled_net','params','colors');
save('temp_gswitch_data','net','rescaled_net','params','col');

function dX = scaled_switchgrad(x)
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
mu = x(1);
nu = x(2);
kappa = x(3);

m = mu;
n = nu*80;
d = kappa*4;

%
dmu = (a0*gamma_0 + a*k0*d)/(gamma_0 + k0*d) - gamma_m*m;
dnu = b*m - gamma_n*n - 2*k1*n^2 + 2*gamma_1*d;
dkappa = k1*n^2 - gamma_1*d;
dX = [dmu;dnu/80;dkappa/4];

end


%function dX = switchgrad(x)
%%%%Parameters
%k0 = 1; %DNA activation
%k1 = 0.0002; %rate of dimerization
%gamma_1 = 2;  %rate of de-dimizeration
%gamma_0 = 50; %DNA inactivation
%gamma_m = 10; %mRNA decay
%gamma_n = 1; %protein decay
%a = 400;	%transcription
%a0 = 0.4;
%b = 40;		%translation
%
%
%%x is a vector in R^3
%m = x(1);
%n = x(2);
%d = x(3);
%
%%
%dm = (a0*gamma_0 + a*k0*d)/(gamma_0 + k0*d) - gamma_m*m;
%dn = b*m - gamma_n*n - 2*k1*n^2 + 2*gamma_1*d;
%dd = k1*n^2 - gamma_1*d;
%
%dX = [dm;dn;dd];
%
%end

function is_far = far(x,net,delta,rho)
%%%NOTE: if input x is part of net, this will return false
%determines if point 'x' is at least delta-far from set of points 'net'
%inputs: x is vector in R^D, net is D X N matrix, N vectors in net
%output: boolean, true if point is far from set

is_far = true;
N = size(net,2);
k = 1;

%Check is x is far away from net
while (is_far && k <= N)
    if rho(x,net(:,k)) < delta
        is_far = false; %break!
    else
        k = k+1;
    end
end

end

function dist = ellipse_dist3(x,y,a,b,c)
%returns distance between x,y on ellipse with axis lengths a,b,c

z = x-y;
dist = sqrt((z(1)/a)^2 + (z(2)/b)^2 + (z(3)/c)^2);

end
