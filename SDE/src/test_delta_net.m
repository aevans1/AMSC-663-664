close all;

%2D example 
x_init = [-1:0.1:2.5];
y_init = [-1:0.1:2];
[X,Y] = meshgrid(x_init,y_init);
init = [X(:)' ; Y(:)'];

rho = @(p_1,p_2) norm(p_1 - p_2);
delta = 0.5;

net = delta_net(init,delta,rho);

scatter(net(1,:),net(2,:),'filled');
%hold on;
