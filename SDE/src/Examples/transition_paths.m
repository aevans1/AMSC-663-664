function transition_paths(atlas_data,transition_params)
%atlas_data is struct with fields
%new_S - new_S from construction.m output, struct of atlas simulator
%L - array of landmarks for atlas corresponding to new_S
%params - setup parameters for atlas, genereated in atlas_driver.m
%('current_atlas.mat' has fields new_S,L from most recent construction.m
%call
%('current_driver.mat' has params struct from most recent current_driver
%call)

%TODO: comments
%general tool for calculating transition paths from atlas output and comparing
%with the original simulator

%NOTE: user can input atlas data and params for this function,
%	or can run it
%with no parameters and have last run atlas .mat files loaded instead

%%%Allow parameter inputs to construcion or load up saved inputs
if ~exist('atlas_data','var')
    try
        fprintf("No input parameters for atlas_data, attempting to load saved params in current_atlas.mat and current_driver.mat \n");
        atlas_data = load('current_atlas.mat');
        data = load('current_driver.mat');
    catch
        fprintf("Error: No saved atlas under current_atlas.mat or current_driver.mat, please generate atlas in construction.m or atlas_driver.m \n");
        return;
    end
end
new_S = atlas_data.new_S;
L = atlas_data.L;
params = data.params;

%%%Check that param struct has non-empty fields
fields = {'S','delta','rho','m','p','t0','d','net_info'};
for i = 1 : length(fields)
    if isempty(getfield(params,fields{i}))
        fprintf("Error: params.%s is empty, please re-initialize param struct with %s value \n", fields{i},fields{i});
        return;
    end
end

%%%Read in delta-net parameters
net_info = params.net_info;
net = net_info.net;
neighbors = net_info.neighbors;
edges = net_info.edges;
deg = net_info.deg;
max_deg = net_info.max_deg;

%%%Read in atlas parameters
delta = params.delta;
d = params.d;
dt = params.dt;
rho = params.rho;
m = params.m;
f = params.f;
dt_sim = params.dt_sim;
S = params.S;


%%%Plot delta-net

%find ambient space dimension
%each column of net is a data point in R^D
D = size(net,1);
if D == 2 || D == 3
    plot_delta_net(net,edges,D);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Define transition regions

%%%Allow parameter inputs to construcion or load up saved inputs
if ~exist('transition_params','var')
    try
        fprintf("No input parameters for transition_params, attempting to load saved params in transition_params.mat \n");
        transition_params = load('transition_params.mat');
    catch
        fprintf("Error: No saved transition parameters under transition_params.mat \n");
        return;
    end
end

%NOTE: including Xzero as a user defined param
%in general, user may want to do several different Xzeros, use a random one,
%but for now assuming that one initial Xzero is fine

%TEMPORARY
regions = transition_params.regions;
dist = transition_params.dist;
dist_sq = dist.^2;
Xzero = transition_params.Xzero;
T = transition_params.T;

%TODO: for more general simulations, use run_opts to specify, e.g: run for
%certain number of transitions rather than certain time
%run_opts = transition_params.run_opts;

%%Assign net points to regions, if within prescribed distance
region_net = zeros(size(net,2),1);
for n = 1 : size(net,2)
    [dmin,imin] = min(sum((net(:,n) - regions).^2,1));
    if dmin < dist_sq
        region_net(n) = imin;
    end
end

num_steps = floor(T/dt_sim);

%%Run simulations



%% Run for original simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Original Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Simulating original paths\n");
%%Split up time interval if time too large to store trajectories in one array
%%Then run code to calculate transitions
original_switches = [];
data_size = num_steps*D; %size of trajectory array for 'num_steps in R^D

if data_size < 1e8
    %%%%simulate one long trajectory [0,T] and calculate transitions
    [~,~,X] = S(Xzero,1,T);
    original_switches = switch_data(X,regions,dist,dt_sim,T);
else
    %%%%chunk up the time interval and calculate switches on the chunks
    %%%%time intervals [increment | increment | ... | rest up time up to T]
    increment = floor(data_size/1e8);
    times = dt_sim*[1e8*ones(1,increment)];
    if mod(data_size,1e8) > 0
        times(end + 1) = dt_sim*mod(data_size,1e8);
    end
    Xinit = Xzero;
    
    %%%%Compute switches up to end of time chunk, then set init at end of chunk
    %and restart, up until time T is reached
    for k = 1:length(times)
        [~,~,X] = S(Xinit,1,times(k));
        switches = 	switch_data(X,regions,dist,dt_sim,times(k));
        original_switches = [original_switches;switches];
        Xinit = X(:,end);
    end
    
end

fprintf("original paths are simulated! \n");
save('original_tswitch.mat','original_switches');




%% Run for learned simulator
num_steps = floor(T/dt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Learned Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Split up time interval if time too large to store trajectories in one array
%%%Then run code to calculate transitions
learned_switches = [];
data_size = num_steps*D; %size of trajectory array for 'num_steps in R^D

%%%Find the initial chart to start reduced simulator
[~,init_chart] = min(sum((Xzero - net).^2,1));
num_nbr = deg(init_chart);
init_L = L(:,1:(m+1)*(num_nbr + 1),init_chart);

%%%Embed initial condition for original simulator,center,then run atlas simulator
[local_L,Yzero] = LMDS(init_L,Xzero,rho,d);
Yzero = Yzero - local_L(:,1);
%%
fprintf("Simulating atlas paths\n");

if data_size < 1e8
    %%%%simulate one long trajectory [0,T] and calculate transitions
    [~,~,~,charts] = learned_simulator(Yzero,1,dt,T,new_S,neighbors,deg,d,delta,net,init_chart);
    %%%Bin the chart list into respective regions, assign 0 if no region
    num_steps = length(charts);
    paths = zeros(1,num_steps);
    for n = 1:num_steps
        paths(n) = region_net(charts(n));
    end
    
    learned_switches = calculate_switches(paths,dt,num_steps);
    
else
    %%%%chunk up the time interval and calculate switches on the chunks
    %%%%time intervals [increment | increment | ... | rest up time up to T]
    increment = floor(data_size/1e8);
    times = dt_sim*[1e8*ones(1,increment)];
    if mod(data_size,1e8) > 0
        times(end + 1) = dt_sim*mod(data_size,1e8);
    end
    
    Yinit = Yzero;
    
    %%%%Compute switches up to end of time chunk, then set init at end of chunk
    %and restart, up until time T is reached
    for k = 1:length(times)
        
        %%%%simulate one long trajectory [0,T] and calculate transitions
        [Yend,~,~,charts] = learned_simulator(Yinit,1,dt,T,new_S,neighbors,deg,d,delta,net,init_chart);
        
        %%%Bin the chart list into respective regions, assign 0 if no region
        num_steps = length(charts);
        paths = zeros(1,num_steps);
        for n = 1:num_steps
            paths(n) = region_net(charts(n));
        end
        
        switches = calculate_switches(paths,dt,num_steps);
        learned_switches = [learned_switches;switches];
        Yinit = Yend;
        init_chart = charts(end);
    end
    
end

fprintf("atlas paths are simulated! \n");
save('learned_tswitch.mat','learned_switches');

%%Save switch data
end

function [switch_times,paths] = switch_data(X,regions,dist,dt,T)
%TODO:comment
%Given trajectory data, calculates paths among regions list,
%finds switch_times betwen regions

num_steps = floor(T/dt);
dist_sq = dist^2;
path_regions = zeros(1,num_steps);
for n = 1:num_steps
    [dmin,imin] = min(sum((regions - X(:,n)).^2,1));
    if dmin < dist_sq
        path_regions(n) = imin;
    else
        path_regions(n) = 0;
    end
end

switch_times = calculate_switches(path_regions,dt,num_steps);
end

function switch_times = calculate_switches(path_regions,dt,num_steps)
%TODO:comment, use transition data file

switch_times=[];
switch_count = 0;
index = path_regions(1);
flag = sign(index);
t = 0;

for n = 2:num_steps
    new_index = path_regions(n);
    if new_index > 0
        if flag == 0
            index = new_index;
            t = 0;
            flag = 1;
        elseif index ~= new_index
            switch_count = switch_count + 1;
            %fprintf("Transition %d \n",switch_count);
            %we've observed a transition!
            %switch_count = switch_count+1
            switch_times = [switch_times;index,new_index,t];
            t = 0;
            index = new_index;
        end
    end
    t = t+dt;
end
end

function plot_delta_net(net,edges,D)
%TODO:comments
%plot a 2d or 3d delta-net

figure; clf; hold on;
set(gca,'Fontsize',20);

if D == 2
    %plot a 2dim delta net, edge by edge
    xlabel('x','Fontsize',20);
    ylabel('y','Fontsize',20);
    for i = 1 : length(edges)
        plot([net(1,edges(i,1)),net(1,edges(i,2))],[net(2,edges(i,1)),net(2,edges(i,2))],'color','k','LineWidth',2);
    end
    drawnow;
    
elseif D == 3
    %plot a 3dim delta net, edge by edge
    xlabel('x','Fontsize',20);
    ylabel('y','Fontsize',20);
    zlabel('z','Fontsize',20);
    for i = 1 : length(edges)
        plot3([net(1,edges(i,1)),net(1,edges(i,2))],[net(2,edges(i,1)),net(2,edges(i,2))] ...
            ,[net(3,edges(i,1)),net(3,edges(i,2))], 'color','k','LineWidth',0.5);
    end
    drawnow;
    view(3);
else
    %can't plot above 3dim or below 2dim!
    fprintf("Error: can only plot 2d or 3d delta-nets \n");
end

end
