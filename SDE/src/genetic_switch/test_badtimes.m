function test_badtimes()

%test bad times

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


sim_time = 10^4;
num_steps = floor(sim_time/dt);

% fprintf("TESTING: setting dt to a smaller value \n");
dt = 0.05;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reduced(ATLAS) Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TESTING: getting running averages of norms of simulations at different
%%%charts
N = size(net,2);

%read in Struct new_S
c = new_S.c;
b = new_S.b;
sigma = new_S.sigma;
T = new_S.T;
mu = new_S.mu;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TESTING: bad steps
%if euler step goes too far ('bad'), collect point and add
%to time spent in 'bad' steps
%bad_points array: badpoint(:,:,n) is array of bad points from chart n,
%                  each column is a differnt data point in R^d

%bad_points =  struct("; %initialize set of bad steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bad_time = 0;
total_time = 0;
bad_time_charts = zeros(1,N);

num_trials = 10;

for trial = 1:num_trials
    
    %pick random chart
    rand_chart = randi(N)
    Xzero = net(:,rand_chart)
    
    fprintf("simulating atlas paths\n");
    
    %%%Find the initial chart to start reduced simulator
    [~,init_chart] = min(sum((Xzero - net).^2,1));
    num_nbr = deg(init_chart);
    init_L = L(:,1:(m+1)*(num_nbr + 1),init_chart);
    
    %%%Embed initial condition for original simulator,center,then run atlas simulator
    [local_L,Yzero] = LMDS(init_L,Xzero,rho,d);
    Yzero = Yzero - local_L(:,1);
    
    y = Yzero;
    i = init_chart; %initial chart

    %simulate until 'switch_count' number of transitions are observed
    for step = 1:num_steps
        
        
        %%%%Code from 'noisy_learned_simulator_step%%%%%%%%%%%%%%%%%%%%%%%%
        %TESTING: explicitly writing out the simulator step, the funciton
        %uses a lot of memory by reloading data structures every step
        %[y,new_chart] = noisy_learned_simulator_step(y,chart,new_S,neighbors,deg,d,dt,delta,noise);
        %[y,new_chart] = learned_simulator_step(y,chart,new_S,neighbors,deg,d,dt,delta);
        
        
        num_nbr = deg(i);
        net_nbr = [i neighbors(i,1:num_nbr)];
        
        
        %%%Find closest chart to y in chart i coords
        C = squeeze(c(:,i,net_nbr)); %collect centers for chart i
        [~,imin] = min(sum((y - C).^2,1));
        j = net_nbr(imin);
        
        if j ~= i
            %TODO: remove squeezes, not necessary for this step
            %%%Set new chart to j
            mu_ij = squeeze(mu(:,i,j));
            mu_ji = squeeze(mu(:,j,i));
            T_ij = squeeze(T(:,:,i,j));
            
            y = T_ij*(y - mu_ij) + mu_ji;
        end
        
        %%%Forward Euler step
        eta = sqrt(dt)*randn(d,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TESTING: no added noise in learned simulator, it's on a
        %different scale than the original simulator
        %y = y + b(:,j)*dt + sigma(:,:,j)*sqrt(noise)*eta;
        y = y + b(:,j)*dt + sigma(:,:,j)*eta;
        
        
        %%%prevent escape from local chart by applying wall function
        if norm(y) > (3/2)*delta
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%TESTING: bad steps
            %if euler step goes too far ('bad'), collect point and add
            %to time spent in 'bad' steps
            bad_time = bad_time + dt;
            bad_time_charts(j) = bad_time_charts(j) + dt;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            y = (1/norm(y))*y*(2*delta - (1/2)*delta*exp(3 - (2/delta)*norm(y)));
            
        end
        i = j;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        total_time = total_time+dt;
        
        
        %NOTE: commented out currently to save on run-time
        %update charts visted list
        %             charts = unique([charts;chart]);
        
        if mod(step,100000) == 0
            fprintf("step %d \n",step);
            bad_time_ratio = bad_time/total_time;
            fprintf("bad_time_ratio = %f \n",bad_time_ratio);
            save('bad_time_data.mat','bad_time_ratio','bad_time_charts');
            toc;
        end
        
        
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TESTING: bad point
    bad_time_ratio = bad_time/total_time;
    fprintf("(bad) time spent stepping to wall: %f \n",bad_time);
    fprintf("time spent: %f \n",total_time);
    fprintf("ratio of bad_time/atlas_time = %f \n",bad_time_ratio);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    save('bad_time_data.mat','bad_time_ratio','bad_time_charts');
    
    %%%%NOTE: commenting out to save running time
    %%%Plots all visited charts by atlas
    %      hold on;
    %      for k = 1:length(charts)
    %          plot3(net(1,charts(k)),net(2,charts(k)),net(3,charts(k)),'.','color','y','Markersize',20);
    %      end
    %      view(3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end
end