load('gswitch_delta_net.mat');
load('atlas_diagnostics.mat');
%%

%%%Plot colors of net_means, more red at point is higher mean of
%%%trajectories of length t_0, more green is lower

color_vals = (t0_list - min(t0_list))./(max(t0_list) - min(t0_list));


%color_vals = (mean_net_distances - min(mean_net_distances))./(max(mean_net_distances) - min(mean_net_distances));


figure; grid; hold on;
for n = 1:size(net,2)
    plot3(net(1,n),net(2,n),net(3,n),'.','Markersize',10,'Color',[color_vals(n),1 - color_vals(n),0]);
    %plot3(net(1,n),net(2,n),net(3,n),'.','Markersize',10,'Color',[C(n,:)]);

end
view(3)
%%

load('bad_time_data.mat');
bad_time_charts = bad_time_charts*1000;
%color_vals = (bad_time_charts - min(bad_time_charts))./(max(bad_time_charts) - min(bad_time_charts))
color_vals = bad_time_charts./max(bad_time_charts)
figure; grid; hold on;
for n = 1:size(net,2)
    plot3(net(1,n),net(2,n),net(3,n),'.','Markersize',10,'Color',[color_vals(n),1 - color_vals(n),0]);
end
view(3)