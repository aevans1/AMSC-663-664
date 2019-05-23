function run_stats_gswitch()
%TODO: comments

close all;

data = load('gswitch_atlas.mat');
new_S = data.new_S;
L = data.L;

data = load('gswitch_driver.mat');
params = data.params;

net_info = params.net_info;
net = net_info.net;
neighbors = net_info.neighbors;
edges = net_info.edges;
deg = net_info.deg;
max_deg = net_info.max_deg;


delta = params.delta;
d = params.d;
dt = params.dt;
rho = params.rho;
m = params.m;
f = params.f;
dt_sim = params.dt_sim;
S = params.S;

%%












%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%T = 10^7; %length of time interval for simulation
%T = 10^5;

dist = 300; %radius of the basin of a equilibrium point
dist_sq = dist.^2;

noise = 100; %optional: noise addedto simulators
max_switches = 1000; %set max number of transitions to track

%%%Define equilibria points and regions of interest for path tracking
%xi - asympt. stable equilibria representing inactive switch
%xa - asympt. stable equilibria representing active switch
xi = [0.0402067142317042; 1.60826856926817; 0.000258652779089588];
xa = [29.3768600805981; 1175.07440322392; 138.079985311206];

global p1
global p2
p1 = xi;
p2 = xa;

regions = [p1,p2];

%%%Assign net points to closest region, if within presciribed distance
region_net = zeros(size(net,2),1);
for n = 1 : size(net,2)
    [dmin,imin] = min(sum((net(:,n) - regions).^2,1));
    if dmin < dist_sq
        region_net(n) = imin;
    end
end

%%%Assign nitial simulation point
Xzero = net(:,1);

%%%Plot Delta Net
% figure;hold on;grid;view(3);
% xlabel('x','Fontsize',20);
% ylabel('y','Fontsize',20);
% set(gca,'Fontsize',20);
% for i = 1 : length(edges)
%     xedge = [net(1,edges(i,1)),net(1,edges(i,2))];
%     yedge = [net(2,edges(i,1)),net(2,edges(i,2))];
%     zedge = [net(3,edges(i,1)),net(3,edges(i,2))];
%     plot3(xedge,yedge,zedge,'color','k','LineWidth',0.25);
% end
% drawnow;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reduced(ATLAS) Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TESTING: getting running averages of norms of simulations at different
%%%charts
N = size(net,2);

%col 1 of running means: running mean of trajectories
%col 2 number of trajectories landed in chart at that index
run_means = zeros(N,1);
run_vars = zeros(N,1);
run_count = zeros(N,1);

fprintf("TESTING learned simulator chart jump norms \n");
noises = [75,60,50,40];
for Q = 1:4
    
    noise = noises(Q);
    fprintf("noise %d \n",noise);
    
    
    
    fprintf("simulating atlas paths\n");
    
    %%%Find the initial chart to start reduced simulator
    [~,init_chart] = min(sum((Xzero - net).^2,1));
    num_nbr = deg(init_chart);
    init_L = L(:,1:(m+1)*(num_nbr + 1),init_chart);
    
    %%%Embed initial condition for original simulator,center,then run atlas simulator
    [local_L,Yzero] = LMDS(init_L,Xzero,rho,d);
    Yzero = Yzero - local_L(:,1);
        
    %%%Compute SDE
    
    %initialize region list, transition tracking
    init_region = region_net(init_chart);
    learned_switches = []; %list of all switches between regions and corresponding times
    charts = []; %list of all charts visited, for plotting
    switch_count = 0; %tracks number of transitions
    tic;
    
    
    %read in Struct new_S
    c = new_S.c;
    b = new_S.b;
    sigma = new_S.sigma;
    T = new_S.T;
    mu = new_S.mu;
    
    %simulate until 'switch_count' number of transitions are observed
    while switch_count < max_switches
        t = 0;
        step = 0; %number of steps iterated
        flag = sign(init_region);
        region = init_region;
        chart = init_chart; %initial chart
        y = Yzero; %initial simulator value in reduced space
        
        
        while region ~= 2
            step = step + 1;
            
            
            %%%%Code from 'noisy_learned_simulator_step
            i = chart;
            
            num_nbr = deg(i);
            net_nbr = [i neighbors(i,1:num_nbr)];
            
            
            %%%Find closest chart to y in chart i coords
            C = squeeze(c(:,i,net_nbr)); %collect centers for chart i
            [~,min_dist] = min(sum((y - C).^2,1));
            j = net_nbr(min_dist);
            
            if j ~= i
                %%%Set new chart to j
                mu_ij = squeeze(mu(:,i,j));
                mu_ji = squeeze(mu(:,j,i));
                T_ij = squeeze(T(:,:,i,j));
                
                y = T_ij*(y - mu_ij) + mu_ji;
            end
            
            %%%Forward Euler step
            eta = sqrt(dt)*randn(d,1);
            y = y + b(:,j)*dt + sigma(:,:,j)*sqrt(noise)*eta;
            
            %%%TESTING: running means, running var%%%%%%%%%%%%%%%%%%%5
            run_count(j) = run_count(j) + 1;
            visits = run_count(j);
            
            run_means(j) = (1/ (visits))*(run_means(j)*(visits -1) + norm(y)); %update running mean
            
            if visits ==1
                run_vars(j) = 0;
            elseif visits == 2
                run_vars(j) = (norm(y) - run_means(j)).^2;
            else
                run_vars(j) = (1/(visits -1))*(run_vars(j)*(visits-2) + (norm(y) - run_means(j)).^2);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%prevent escape from local chart by applying wall function
            if norm(y) > (3/2)*delta
                y = (1/norm(y))*y*(2*delta - (1/2)*delta*exp(3 - (2/delta)*norm(y)));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            new_chart = j;
            
            
            %[y,new_chart] = noisy_learned_simulator_step(y,chart,new_S,neighbors,deg,d,dt,delta,noise);
            %[y,new_chart] = learned_simulator_step(y,chart,new_S,neighbors,deg,d,dt,delta);
            
            chart = new_chart;
            
            %calculate transition path
            new_region = region_net(chart);
            if new_region > 0
                if flag == 0
                    region = new_region;
                    t = 0;
                    flag = 1;
                else
                    if region ~= new_region
                        
                        %we've observed a transition!
                        switch_count = switch_count + 1;
                        learned_switches = [learned_switches;region,new_region,t];
                        t = 0;
                        region = new_region;
                        
                        %TESTING
                        filename = ['running_stats','_','noise',num2str(noise),'.mat'];
                        save(filename,'run_means','run_vars','run_count');
                        
                        fprintf("transition %d \n",switch_count);
                        toc;
                    end
                end
            end
            t = t+dt;
            
            
            %update charts visted list
            %         charts = unique([charts;chart]);
            
            if mod(step,1000000) == 0
                fprintf("step %d \n",step);
                % 		fprintf("step %d of %f \n",step,N);
                toc;
            end
        end
    end
    atlas_time = toc;
    
    
    
    filename = ['running_stats','_','noise',num2str(noise),'.mat'];
    save(filename,'run_means','run_vars','run_count');
    
end
%TESTING: noise values for chart norm jumps

fprintf("atlas paths are simulated! \n");

%Save info relevant for plotting mean switching times,
%and save rest of information for archival purposes

end



function dX = switchgrad(x)
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