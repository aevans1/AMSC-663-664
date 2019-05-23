data = load('gswitch_data.mat');

delta = 2.0;
init = data.net;
rho = @(x,y) norm(x - y);
delta_net(init,delta,rho);
