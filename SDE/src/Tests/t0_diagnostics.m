load('current_delta_net.mat');
load('current_driver.mat');
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

D = size(net,1);

%% Create Graph Structure and Order Nodes via graph search

my_node = 1;  %Starting node for the graph search

%varying t0 values will be assigned according to a neighbor search
%in this case a breadth first search through the delta-net
delta_graph = graph(edges(:,1),edges(:,2));
bf_ordering = bfsearch(delta_graph,my_node);

%% Assign new t0 values for net points

% -Traverse the graph by some search algorithm (bfS/bfS)
% 	-Given a point x,
% -set t_0 to average of t_0 for neighbors
% 			-if no neighbors, set to some default
% 		-if delta_avg > 1.1 delta
% 			-t_0:= t_0 / 2
% 		-if delta_avg < 0.2 delta
% 			-t_0:=2*delta 

mean_net_distances = zeros(1,length(bf_ordering));
t0_list = zeros(1,length(bf_ordering));
for k = 1:length(bf_ordering)
    fprintf("Step %d \n", k);
    n = bf_ordering(k);
    
    
    %Set t0 value average of neighbors mean_net_dist
    num_nbr = deg(n);
    nbr = neighbors(n,1:num_nbr);
    nbr_t0 = t0_list([nbr]);
    nbr_t0 = nbr_t0(nbr_t0 > 0); %only use neighbors visited so far
    
    %set t0 value at node to mean of neighbors if neighbors have been
    %visited, else set at default t0
    if ~isempty(nbr_t0)
        new_t0 = mean(nbr_t0); %avg over neighbors
    else
        new_t0 = t0; %initialize
    end
        
    %%%Simulate p paths around net point y_n
    X(:,1:p,n) = S(net(:,n),p,new_t0);
    
    mean_net_dist = test_mean_net_dist(n,D,net,delta,new_t0,X);
    
    %%%Varying t_o values based on how far points go in dynamics
    if mean_net_dist >1.5*delta
        %NOTE:  > 1.5 delta means wall function will be applied ALOT, not good
        
        new_t0 = new_t0/2;
        X(:,1:p,n) = S(net(:,n),p,new_t0);
        mean_net_dist = test_mean_net_dist(n,D,net,delta,new_t0,X);
    elseif mean_net_dist < 0.5*delta
        new_t0 = 2*new_t0;
        X(:,1:p,n) = S(net(:,n),p,new_t0);
        mean_net_dist = test_mean_net_dist(n,D,net,delta,new_t0,X);
    end
    
    t0_list(k) = new_t0;
    mean_net_distances(k) = mean_net_dist;
end
    
save('atlas_diagnostics.mat','mean_net_distances','t0_list');

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

