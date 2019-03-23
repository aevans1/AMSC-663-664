function MakeDeltaNet(delta,pot,grad)
N = 100;
xmin = -1.0;
xmax = 2.5;
ymin = -1.0;
ymax = 2;
x = linspace(xmin,xmax,N);
y = linspace(ymin,ymax,N);
[xg,yg] = meshgrid(x,y); 
for i = 1 : N
    for j = 1 : N
        U(i,j) = pot([xg(i,j);yg(i,j)]);
    end
end
% plot the level sets of the potential
figure(1);
clf; hold on;
contour(x,y,U,linspace(0,3.7,20),'Linewidth',1);
% run a long trajectory
dt = 0.005;
nsteps = 800000;
dw = sqrt(dt)*randn(2,nsteps);
x = 0.5*[xmin + xmax;ymin + ymax];
net(:,1) = x;
nnet = 1;
d2 = delta^2;
for i = 1 : nsteps
    x = x + grad(x)*dt + 10*dw(:,i);
    d = sum((net - x*ones(1,nnet)).^2,1);
    if min(d) > d2 & pot(x) < 10
        nnet = nnet + 1;
        net(:,nnet) = x;
        plot(x(1),x(2),'o','color','k');
        drawnow;
    end
end
%% make the delta-net into a network and remove isolated vertices
del2 = 2*delta;
E = []; % edges of the delta-net
for i = 1 : nnet - 1
    for j = i + 1 : nnet
        if norm(net(:,i) - net(:,j)) < del2
            plot(net(1,[i,j]),net(2,[i,j]),'color','k','Linewidth',1);
            E = [E;[i,j]];
        end
    end
end

% degrees and maximal degree
n0 = 0;
Nnet = size(net,2);
deg = zeros(Nnet,1);
isolated = [];
for k = 1 : Nnet
    ind = find(E(:,1) == k | E(:,2) == k);
    if isempty(ind)
        n0 = n0 + 1;
        isolated(n0) = k;
    else
        deg(k) = length(ind);
    end
end
% remove isolated vertices
if ~isempty(isolated)
    net(:,isolated) = [];
    deg(isolated) = [];
    Nnet = Nnet - n0;
end
maxdeg = max(deg); % the maximal degree

E = []; % edges of the delta-net
for i = 1 : nnet - 1
    for j = i + 1 : nnet
        if norm(net(:,i) - net(:,j)) < del2
            plot(net(1,[i,j]),net(2,[i,j]),'color','k','Linewidth',1);
            E = [E;[i,j]];
        end
    end
end

maxdeg = max(deg); % the maximal degree
% find neighbors
% vertices k and j are neighbors, k ~ j, iff |net(:,k)-net(:,j)| < 2*delta
nei = zeros(Nnet,maxdeg);
for k = 1 : Nnet
    ind = find(E(:,1) == k);
    if ~isempty(ind)
        l1 = length(ind);
        nei(k,1 : l1) = E(ind,2)';
    else
        l1 = 0;
    end
    ind = find(E(:,2) == k);
    if ~isempty(ind)
        l2 = length(ind);
        nei(k,l1 + 1 : l1 + l2) = E(ind,1)';
    end
end




fsize = 20;
colorbar
daspect([1,1,1]);
xlabel('x_1','Fontsize',fsize);
ylabel('x_2','Fontsize',fsize);
set(gca,'Fontsize',fsize);

save('delta_net_CM2D.mat','net','E','deg','maxdeg','nei','xg','yg','U');
end



