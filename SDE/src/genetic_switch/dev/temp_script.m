%run_atlas_tpaths_gswitch;
%makestats_tpaths_gswitch;

%%
data = load('gswitch_atlas.mat');
new_S = data.new_S;
L = data.L;

data = load('gswitch_driver.mat');
params = data.params;

net_info = params.net_info;
net = net_info.net;
edges = net_info.edges;

delta = params.delta;
d = params.d;

%%%Define equilibria points and regions of interest for path tracking

%xi - asympt. stable equilibria representing inactive switch
%xa - asympt. stable equilibria representing active switch
%xs - saddle
xi = [0.0402067142317042; 1.60826856926817; 0.000258652779089588];
xa = [29.3768600805981; 1175.07440322392; 138.079985311206];
xs = [10.5829; 423.3173; 17.9198];

global p1
global p2
p1 = xi;
p2 = xa;

regions = [p1,p2];


dist = 200;
dist_sq = dist.^2;

%%%Assign net points to closest region, if within presciribed distance
region_net = zeros(size(net,2),1);
for n = 1 : size(net,2)
    [dmin,imin] = min(sum((net(:,n) - regions).^2,1));
    if dmin < dist_sq
        region_net(n) = imin;
    end
end

figure;hold on;grid;view(3);
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
set(gca,'Fontsize',20);
% for i = 1 : length(edges)
%     xedge s= [net(1,edges(i,1)),net(1,edges(i,2))];
%     yedge = [net(2,edges(i,1)),net(2,edges(i,2))];
%     zedge = [net(3,edges(i,1)),net(3,edges(i,2))];
%     plot3(xedge,yedge,zedge,'color','k','LineWidth',0.25);
% end
% drawnow;
for n = 1:size(net,2)
    switch region_net(n)
        case 0 
            color_val = 'k';
        case 1 
            color_val = 'b';
        case 2 
            color_val = 'g';        
    end
    plot3(net(1,n),net(2,n),net(3,n),'.','color',color_val);
end


%plot equilibria and saddle
plot3(xi(1),xi(2),xi(3),'.','color','g','Markersize',20);
plot3(xa(1),xa(2),xa(3),'.','color','g','Markersize',20);
plot3(xs(1),xs(2),xs(3),'.','color','g','Markersize',20);



%%

num_points = 500;
box_init = [rand(1,num_points)*40;rand(1,num_points)*1500;rand(1,num_points)*200];
box_new = zeros(3,num_points);
noise = 1;
%adjusting dt_sim
dt_sim = 0.01;
N = 10^2;
num_steps = floor(N/dt_sim);

for i = 1:num_points
    x = box_init(:,i);
    plot3(x(1),x(2),x(3),'.','color','r','LineWidth',5.0);
    
    for n = 1:num_steps
        dW = sqrt(dt_sim)*randn(3,1); %Brownian random increment in R^D
        x = x + switchgrad(x)*dt_sim + sqrt(noise)*dW;
    end
    plot3(x(1),x(2),x(3),'.','color','m','LineWidth',5.0);
    box_new(:,i) = x;
end



%%

% 
% noise = 100;
% 
% %adjusting dt_sim
% dt_sim = 0.01;
% N = 10^6;
% num_steps = floor(N/dt_sim);
% x = Xzero;
% 
% for n = 1:num_steps
%     dW = sqrt(dt_sim)*randn(3,1); %Brownian random increment in R^D
%     
%     %TESTING: for plotting blow-ups
%     increment = switchgrad(x)*dt_sim + sqrt(noise)*dW;
%     x = x + increment;
%     
%     if mod(n,100) == 0
%         %      if sum(abs(increment.^2),2) > 5
%         plot3(x(1),x(2),x(3),'.','color','m','LineWidth',5.0);
%         %         fprintf("plotting!\n")
%         %          drawnow;
%     end
%     
%     
%     if mod(n,100000000) == 0
%         fprintf("step %d \n",n);
%         % 		fprintf("step %d of %f \n",step,N);
%         toc;
%     end
%     
% end
% 
% 






%%
% for n = 1:size(net,2)
%
%     point = net(:,n);
%     region_val = region_net(n);
%     switch region_val
%         case 0
%         case 1
%             plot3(point(1),point(2),point(3),'.','color','r','Markersize',20);
%         case 2
%             plot3(point(1),point(2),point(3),'.','color','b','Markersize',20);
%
%     end
%
% end

%%
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
