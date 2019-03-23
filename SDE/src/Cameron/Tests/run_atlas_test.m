function run_atlas_test()
% runs atlas on the original 2D space, hence avoids the use of
% interchart switch maps and estimates for b and sigma 
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
fprintf('delta = %d, Nnet = %d\n',delta,Nnet);
%%
% obtain exact b and sigma at the delta-net points
pot = @(x)CMpot2D(x);
fun = @(x)(-CMgrad2D(x)); % dx = -grad(U(x))dt + dW
for k = 1 : Nnet
    b(:,k) = fun(net(:,k));
end
% sigma is taken to be I as all charts    
I = eye(2);


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
well = zeros(Nnet,1);
wrad = 0.25; % radius of the ball surrounding the wells
for k = 1 : Nnet
    [dmin,imin] = min([norm(net(:,k) - p1),norm(net(:,k) - p2),norm(net(:,k) - p3)]);
    plot(net(1,k),net(2,k),'.','color',col(imin),'MarkerSize',30);
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
fprintf('dt = %d\n',dt);
dw = sqrt(dt)*randn(2,Nsteps);
% for the Wall-map 
% W(x) =(x/|x|)*(2*Delta - (delta/2)*exp(3 - (2/delta)*|x|))
d32 = delta*1.5;
d2 = 2*delta;
d05 = 0.5*delta; % delta/2
ed05 = 2/delta; % 2/delta

chart = 1;
x = net(:,chart); % start with the first delta-net point
t = 0;
iw = well(chart);
flag = sign(iw);

for n = 1 : Nsteps
    neibs = [chart,nei(chart,1 : deg(chart))];
    C = net(:,neibs); 
    [xmin,imin] = min(sum((x*ones(1,length(neibs)) - C).^2,1));
    newchart = neibs(imin);
    ch(n) = newchart;
    iwnew = well(newchart);
    if iwnew > 0 % got into some well
        if flag == 0 % this is the first well
            iw = iwnew;
            flag = 1;
            t = 0;
        else
            if iw ~= iwnew
                Tswitch = [Tswitch; [iw,iwnew,t]];
                iw = iwnew;
                t = 0;
            end
        end
    end
    % forward Euler step
    eta = randn(2,1);
    t = t + dt;
    x = x + b(:,newchart)*dt + dw(:,n);
    % prevent escape from local chart
    y = x - net(:,chart);
    normx = norm(y);
    if normx > d32
        y = y*(d2 - d05*exp(3 - ed05*normx))/normx;
        x = y + net(:,chart);
    end
    chart = newchart; 
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



