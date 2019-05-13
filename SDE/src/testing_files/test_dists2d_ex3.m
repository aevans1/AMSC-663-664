% load('gswitch_delta_net.mat');
% load('distance_checks.mat');
load('current_delta_net.mat');
load('current_atlas.mat');
%%

%%%Plot colors of net_means, more red at point is higher mean of
%%%trajectories of length t_0, more green is lower
net_means = atlas_diagnostics.mean_net_distances;

net_means = (net_means > 0.4).*(net_means);

color_vals = (net_means - min(net_means))./(max(net_means) - min(net_means));


figure; grid; hold on;
for n = 1:size(net,2)
    plot(net(1,n),net(2,n),'.','Markersize',10,'Color',[color_vals(n),1 - color_vals(n),0]);
    %plot3(net(1,n),net(2,n),net(3,n),'.','Markersize',10,'Color',[color_vals(n),1 - color_vals(n),0]);
end
%view(3);
%%
%%%Redo above, but with values under 1.5*delta (1.5 * 5 = 7.5) set to 0
net_means = (net_means > 7.5).*(net_means);
color_vals = (net_means - min(net_means))./(max(net_means) - min(net_means));

figure; grid; hold on;
for n = 1:size(net,2)
    plot3(net(1,n),net(2,n),net(3,n),'.','Markersize',10,'Color',[color_vals(n),1 - color_vals(n),0]);
end
view(3);
