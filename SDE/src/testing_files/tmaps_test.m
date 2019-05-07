function tmaps_test()

data = load('current_delta_net.mat');
net = data.net;
neighbors = data.neighbors;
edges = data.edges;
deg = data.deg;
max_deg = data.max_deg;

data = load('current_atlas.mat');
new_S = data.new_S;
T = new_S.T;
mu = new_S.mu;
embed_L = data.embed_L;
L = data.L;

params = data.params;
delta = params.delta;
d = params.d;
dt = params.dt;
rho = params.rho;
m = params.m;
f = params.f;
dt_original = params.dt_original;
S = params.S;



L1new = embed_L(:,:,1);
L22new = embed_L(:,:,22);

figure;hold on; grid;
col1 = 'b';
col2 = 'r';
plot(L1new(1,1:m+1),L1new(2,1:m+1),'.','color',col1,'Markersize',20);
plot(L1new(1,1),L1new(2,1),'.','color',col1,'Markersize',40);
plot(L1new(1,m+2:2*(m+1)),L1new(2,m+2:2*(m+1)),'.','color',col2,'Markersize',20);
plot(L1new(1,m+2),L1new(2,m+2),'.','color',col2,'Markersize',40);
daspect([1,1,1]);
  

i = find([22 neighbors(22,:)] == 1);

figure; hold on; grid;
col1 = 'r';
col2 = 'b';
plot(L22new(1,1:m+1),L22new(2,1:m+1),'.','color',col1,'Markersize',20);
plot(L22new(1,1),L22new(2,1),'.','color',col1,'Markersize',40);
plot(L22new(1,(i-1)*(m+1) + 1:(i)*(m+1)),L22new(2,(i-1)*(m+1) + 1:(i)*(m+1)),'.','color',col2,'Markersize',20);
plot(L22new(1,(i-1)*(m+1) + 1),L22new(2,(i-1)*(m+1) + 1),'.','color',col2,'Markersize',40);
daspect([1,1,1]);
   

ch1pt1 = L1new(:,1:m+1);
ch1pt2 = L1new(:,m+2:2*(m+1));
ch2pt1 = L22new(:,(i-1)*(m+1) + 1)
ch2pt2 = L22new(:,1:m+1);

T12 = squeeze(T(:,:,1,22));
T21 = squeeze(T(:,:,22,1));
mu12 = mu(:,1,22);
mu21 = mu(:,22,1);

map12pt1 = T12*(ch1pt1 - mu12) + mu21;

plot(map12pt1(1,:),map12pt1(2,:),'.','color','k','Markersize',30);
 
end
