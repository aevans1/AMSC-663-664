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

%TEMPORARY
D = 3;


net_info= params.net_info;
net = net_info.net;
neighbors = net_info.neighbors;
edges = net_info.edges;
deg = net_info.deg;
max_deg = net_info.max_deg;



%%
%Plot Delta Net
figure(1);hold on;grid;view(3);
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
set(gca,'Fontsize',20);
for i = 1 : length(edges)
    xedge = [net(1,edges(i,1)),net(1,edges(i,2))];
    yedge = [net(2,edges(i,1)),net(2,edges(i,2))];
    zedge = [net(3,edges(i,1)),net(3,edges(i,2))];
    plot3(xedge,yedge,zedge,'color','k','LineWidth',0.05);
end
drawnow;
%%
%highlight a node and neighbors
figure(1);hold on;
my_node = 1155;
num_nbr = deg(my_node);
nbr = neighbors(my_node,1:num_nbr);
for i = 1 : num_nbr
    
    xedge = [net(1,my_node),net(1,nbr(i))];
    yedge = [net(2,my_node),net(2,nbr(i))];
    zedge = [net(3,my_node),net(3,nbr(i))];
    plot3(xedge,yedge,zedge,'color','g','LineWidth',1);
end


plot3([net(1,my_node),net(1,my_node)],[net(2,my_node),net(2,my_node)],[net(3,my_node),net(3,my_node)],'.','color','b','MarkerSize',10);

%%

weights = [];
for n = 1:size(edges,1)
   x = edges(n,1);
   y = edges(n,2);
   weights(n) = norm(net(:,x) - net(:,y));
end

%%
delta_graph = graph(edges(:,1),edges(:,2));

figure
graph_plot = plot(delta_graph);
layout(graph_plot,'force');
%view(3)

%figure out on graph where one 'endpoint' of the graph is to start the
%search

%%
my_node = 2645; %found from analyzing graph_plot from above, 'end node'
%varying t0 values will be assigned according to a neighbor search
%in this case a breadth first search through the delta-net
delta_graph = graph(edges(:,1),edges(:,2));
bf_ordering = bfsearch(delta_graph,my_node);



% -Traverse the graph by some search algorithm (bfS/bfS)
% 	-Given a point x,
% -set t_0 to average of t_0 for neighbors
% 			-if no neighbors, set to some default
% 		-if delta_avg > 1.1 delta
% 			-t_0:= t_0 / 2
% 		-if delta_avg < 0.2 delta
% 			-t_0:=2*delta
% -(learned simulator) traverse graph again
% 	-run learned simulator for certain number of steps (e.g 5), see if wall funciton is applied
% 		-adjust h for this chart so that points don’t hit wall in that prescribed number of steps
% 

mean_net_distances = zeros(1,length(bf_ordering));
t0_list = zeros(1,length(bf_ordering));
for k = 1:length(bf_ordering)
    
    fprintf("Step %d \n", k);
    
    n = bf_ordering(k);
    
    
        %find average of neighbors mean_net_dist
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
    
    mean_net_dist = test_mean_njkket_dist(n,D,net,delta,new_t0,X);
    
    %%%Varying t_o values based on how far points go in dynamics
    if mean_net_dist >1.5*delta
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

function [dist_to_net] = test_mean_net_dist(n,D,net,delta,t0,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO: replace euclidean norm with more general norm
%t0 simulation time should be chosen so that
%E[X - net(:,n)] ~ delta, likewise for landmarks L
if D == 1
    dist_to_net =  mean(abs(X(:,:,n) - net(:,n)),2);
else
    dist_to_net =  mean(vecnorm(X(:,:,n) - net(:,n)),2);
end
fprintf("avg distance from X to net point n:	%f \n",dist_to_net);
fprintf("	delta: 	%f\n",delta);
fprintf("	t_0: 	%f\n",t0);

end
