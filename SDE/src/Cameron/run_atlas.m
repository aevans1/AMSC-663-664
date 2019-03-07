function run_atlas()
global p1 p2 p3
p1 = [0;0];
p2 = [1.5;0];
p3 = [0.8;1.05];
ppp = [p1,p2,p3];
% save('LearnedSimulator_CM2D.mat','T','mu','c','b','sig');
%% OUTPUT of learn_atlas
% The number of charts is Nnet
% (1) For map from chart k to chart j:
% T = zeros(d,d,Nnet,Nnet)are  Nnet*Nnet d-by-d matrices
% mu = zeros(2*m1,Nnet,Nnet) are Nnet*Nnet d-by-1 vectors
% (2) c = zeros(d,Nnet,Nnet); % c_{k,j} = Phi_k(y_j), 
% where Phi_k is the embedding of the k-th chart
% (3) Drifts and covariance matrices for simulator at chart k:
% b = zeros(d,Nnet);
% sig = zeros(d,d,Nnet);
%% Read data
dd = load('LearnedSimulator_CM2D.mat');
b = dd.b;
sig = dd.sig;
T = dd.T;
mu = dd.mu;
c = dd.c;
delta = dd.delta;

data = load('delta_net_CM2D.mat');
net = data.net; % coordinates of vertices of the delta-net
E = data.E; % edges of the delta-net
xgrid = data.xg;
ygrid = data.yg;
Ugrid = data.U;
deg = data.deg;
maxdeg = data.maxdeg;
nei = data.nei;


Nnet = size(net,2);

figure(1); clf; hold on;
contour(xgrid,ygrid,Ugrid,linspace(0,3.7,20),'Linewidth',1);
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
set(gca,'Fontsize',20);
for i = 1 : length(E)
    plot([net(1,E(i,1)),net(1,E(i,2))],[net(2,E(i,1)),net(2,E(i,2))],'color','k','LineWidth',2);
end
drawnow;

%% Determine which points of delta-net relate to which well
col = ['r','g','b'];
fun = @(x)(-CMgrad2D(x)); % dx = -grad(U(x))dt + dW
odefun = @(t,x)fun(x);
Tspan = [0,100];
options = odeset('Events',@events);
basin = zeros(Nnet,1);
well = zeros(Nnet,1);
wrad = 0.25; % radius of the ball surrounding the wells
for k = 1 : Nnet
    [~,Y] = ode45(odefun,Tspan,net(:,k),options);
    yf = Y(end,:)';
    [dmin,imin] = min([norm(yf - p1),norm(yf - p2),norm(yf - p3)]);
    plot(net(1,k),net(2,k),'.','color',col(imin),'MarkerSize',30);
    basin(k) = imin;
    if norm(net(:,k) - ppp(:,imin)) < wrad
        well(k) = imin;
        fprintf('k = %i, well = %i\n',k,well(k));
    end
end
    
    

%% Run ATLAS
Nsteps = 1000000; % the number of atlas steps
Nwell = 3;
Tswitch = []; % switching times between wells
dt = 0.2*delta^2;
sqdt = sqrt(dt);
% for the Wall-map 
% W(x) =(x/|x|)*(2*Delta - (delta/2)*exp(3 - (2/delta)*|x|))
d32 = delta*1.5;
d2 = 2*delta;
d05 = 0.5*delta; % delta/2
ed05 = 2/delta; % 2/delta

chart = 1;
x = squeeze(c(:,chart,chart)); % start with the first delta-net point
tstart = 0;
t = 0;
wellprev = well(chart);
if wellprev > 0 
    wellflag = 1;
else 
    wellflag = 0;
end
ch = zeros(Nsteps,1);
for n = 1 : Nsteps
    neibs = [chart,nei(chart,1 : deg(chart))];
    C = squeeze(c(:,chart,neibs)); 
    [xmin,imin] = min(sum((x*size(C,2) - C).^2,1));
    newchart = neibs(imin);
    ch(n) = newchart;
    if chart ~= newchart
        sh1 = squeeze(mu(:,chart,newchart));
        sh2 = squeeze(mu(:,newchart,chart));
        TT = squeeze(T(:,:,chart,newchart));
        x = ((x - sh1)'*TT)' + sh2;
        if well(newchart) > 0 & well(newchart) ~= wellprev
            w1 = wellprev;
            w2 = well(newchart);
            Tswitch = [Tswitch;w1,w2,t - tstart];
            tstart = t;
            wellprev = w2;
        end
    end
    % forward Euler step
    eta = randn(2,1);
    t = t + dt;
    x = x + b(:,newchart)*dt + squeeze(sig(:,:,newchart))*eta*sqdt;
    % prevent escape from local chart
    normx = norm(x);
    if normx > d32
        x = x*(d2 - d05*exp(3 - ed05*normx));
    end
    chart = newchart;    
    if wellflag == 0
        wellprev = well(chart);
        if wellprev > 0
            wellflag = 1;
        end
    end
%     fprintf('n = %i, chart = %i, basin = %i, well = %i\n',n,chart,basin(chart),well(chart));    
end
save('TswitchCM2D.mat','Tswitch','ch','Nsteps','dt','wrad');
% mark all visited delta net points
plot(net(1,ch),net(2,ch),'.','color','y','Markersize',20);
end
%%
function [position,isterminal,direction] = events(t,x)
global p1 p2 p3
position = max(0,min([norm(x - p1),norm(x - p2),norm(x - p3)]) - 0.1);
isterminal = 1; % terminate integration
direction = 0;
end



