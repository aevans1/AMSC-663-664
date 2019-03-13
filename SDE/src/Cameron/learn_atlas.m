function learn_atlas()
%%
% Crosskey and Maggioni (2017) ATLAS, SIAM MMS code from Fig. 4, page 119
% 2D example with a 3-well smooth potential, page 142
% CM2D, CMpot2D, CMgrad2D
%% OUTPUT
% The number of charts is Nnet
% (1) For map from chart k to chart j:
% T = zeros(d,d,Nnet,Nnet)are  Nnet*Nnet d-by-d matrices
% mu = zeros(2*m1,Nnet,Nnet) are Nnet*Nnet 2*(m + 1)-by-1 vectors
% (2) c = zeros(d,Nnet,Nnet); % c_{k,j} = Phi_k(y_j), 
% where Phi_k is the embedding of the k-th chart
% (3) Drifts and covariance matrices for simulator at chart k:
% b = zeros(d,Nnet);
% sig = zeros(d,d,Nnet);
%%
%% Data and parameters
delta = 0.2;
pot = @(x)CMpot2D(x);
fun = @(x)(-CMgrad2D(x)); % dx = -grad(U(x))dt + dW
 MakeDeltaNet(delta,pot,fun); % created and saved
data = load('delta_net_CM2D.mat');
net = data.net; % coordinates of vertices of the delta-net
E = data.E; % edges of the delta-net
xgrid = data.xg;
ygrid = data.yg;
Ugrid = data.U;
deg = data.deg;
maxdeg = data.maxdeg;
nei = data.nei;
%% graphics
figure(2); clf; hold on;
plot(net(1,:),net(2,:),'.','Markersize',20);
le = length(E);
for i = 1 : le
    plot([net(1,E(i,1)),net(1,E(i,2))],[net(2,E(i,1)),net(2,E(i,2))],'color','k','LineWidth',2);
end
contour(xgrid,ygrid,Ugrid,linspace(0,3.7,20),'Linewidth',1);
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
set(gca,'Fontsize',20);
drawnow;
%%
% parameters
dim = size(net,1); % dimension of the whole space where manifold M is embedded
d = 2; % dimension of the manifold M
Nnet = size(net,2); % the number of points in the delta-net
t0 = delta^2; % I added the factor 0.1; the simulation time forgetting landmarks
m = 10; % m + 1 = the number of landmarks
p = 10000; % 10000 number of paths to determine drift and diffusion
%% create landmarks
m1 = m + 1;
for k = 1 : Nnet
    A(:,1,k) = net(:,k);
    a = simulator(net(:,k),m,t0,fun); % a is d-by-m matrix where x\in R^d
    A(:,2 : m1,k) = a;
end
% A(:,1,k) is the k-th point of the delta-net
% A(:,2 : m1,k) are the landmarks around k-th point of the delta-net

% Plot landmarks and connect them with straight line segments with the
% corresponding points of the delta-net
col = jet(Nnet);
for k = 1 : Nnet
    for j = 1 : m1
%         fprintf('k = %d, j = %d, net = [%d,%d], A = [%d,%d]\n',k,j,net(1,k),net(2,k),A(1,j,k),A(2,j,k));
        if j >= 2
            plot([A(1,1,k),A(1,j,k)],[A(2,1,k),A(2,j,k)],'color',col(k,:));
        else
            Ap = squeeze(A(:,:,k));
            plot(Ap(1,:),Ap(2,:),'*','color','r');
        end
    end
end
%% simulate paths for estimating drift and diffusion for each net point
% note: k ~ j if |net(:,k) - net(:,j)| < 2*delta
del2 = 4*delta^2;
pt0 = p*t0;
c = zeros(d,Nnet,Nnet); % c_{k,j} = Phi_k(y_j), Phi_k is the embedding of the k-th chart
Lnew = zeros(d,(maxdeg + 1)*m1,Nnet); % L(:,:,k) = Phi_k(Lk)
T = zeros(d,d,Nnet,Nnet);
mu = zeros(d,Nnet,Nnet);
b = zeros(d,Nnet);
sig = zeros(d,d,Nnet);
for k = 1 : Nnet
    fprintf('k = %d\n',k);
    x = simulator(net(:,k),p,t0,fun);
    Lk = squeeze(A(:,:,k));
    index = [];
    % index are the indices of net(:,j) in Lk
    for neib = 1 : deg(k)
        % find points for landmarks within radius 2*delta of net(:,k)
        j = nei(k,neib);
        A1 = squeeze(A(:,:,j));
        index = [index;[j,size(Lk,2) + 1]];
        Lk = [Lk,A1]; % Lk = union_{i ~ k} A_i
    end
    % Lk is the union of k-th delta-net point together with its landmarks 
    % and all delta-net points lying with distance 2*delta from the k-th
    % delta-net point + their landmarks
    [Lknew,xnew] = LMDS(Lk,x,d); % d is the desired output dimension
    % find Phi_k(net(:,j)) where Phi_k is the embedding map of the k-th
    % chart
     nneib = size(index,1);
    % shift coordinates so that c(:,k,k) = 0
    shift = Lknew(:,1); % place the image of k-th delta-net point into the origin
    Lknew = Lknew - shift*ones(1,size(Lknew,2));
    xnew = xnew - shift*ones(1,size(xnew,2));
    for neib = 1 : nneib
        c(:,k,index(neib,1)) = Lknew(:,index(neib,2)); 
    end
    c(:,k,k) = zeros(d,1);
%     figure; hold on
%     plot(x(1,:),x(2,:),'.','Markersize',20,'color',[0.5 0 0.5]);
%     plot(Lk(1,:),Lk(2,:),'o','color','m');
%     figure;
%     hold on
%     plot(xnew(1,:) + net(1,k),xnew(2,:) + net(2,k),'.','Markersize',20,'color',[0.5 0 0.5]);
%     plot(Lknew(1,:) + net(1,k),Lknew(2,:) + net(2,k),'o','color','m');
    Lnew(:,1 : m1*(deg(k) + 1),k) = Lknew;
    % find drift and diffusion
    b(:,k) = sum(xnew,2)/pt0;
    aux = cov(xnew')/t0;
    [evec,eval] = eig(aux);
    sig(:,:,k) = evec*sqrt(eval)*evec';
    % compute switching maps
    for neib = 1 : nneib
        j = index(neib,1);
        i = index(neib,2);
        if j < k
            Lkjnew = [Lknew(:,1 : m1),Lknew(:,i : i + m1 - 1)]';
            Ljnew = Lnew(:,1 : m1*deg(j),j);
            ik = find(nei(j,:) == k);
%             size(Lknew)
%             Lknew
            Ljknew = [Ljnew(:,1 : m1),Ljnew(:,ik : ik + m1 - 1)]';
%             fprintf('j = %d, k = %d\n',j,k);
%             size(Lkjnew)
%             size(Ljknew)
%             Lkjnew
%             Ljknew
            % Lkjnew and Ljknew are 2*m1-by-d matrices
            % find their row means
            mukj = mean(Lkjnew,1);
            mujk = mean(Ljknew,1);
            Tkj = (MPpseudo(Lkjnew - ones(size(Lkjnew,1),1)*mukj))*(Ljknew - ones(size(Lkjnew,1),1)*mujk);
            mu(:,k,j) = mukj';
            mu(:,j,k) = mukj';
            T(:,:,k,j) = Tkj;
            T(:,:,j,k) = Tkj;
        end
    end
end
save('LearnedSimulator_CM2D.mat','delta','T','mu','c','b','sig');
end
%% Moore-Penrose pseudoinverse
function B = MPpseudo(A)
[m,n] = size(A);
flag = 0;
if m > n
    r = rank(A);
    if r == n
        B = (A'*A)\(A');
        flag = 1;
    end
end
if flag == 0
    fprintf('Error in MPinverse: columns of A are linearly dependent\n');
    B = 0;
end
end
        


