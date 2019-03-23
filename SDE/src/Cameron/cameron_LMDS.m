function [newlmarks,newdata] = cameron_LMDS(lmarks,data,k)
%% Input
% lmarks = landmarks = d-ny-n matrix
% data = other data points, d-by-Ndata matrix
% k = the desired number of output dimensions;
%% setup
fprintf('In LMDS\n');
Ndata = size(data,2); % data is d-by-N
n = size(lmarks,2); % the number of landmarks
N = Ndata + n; % the size of the whole data set;
% %% plot input data
% figure(3); clf; hold on;
% plot(data(1,:),data(2,:),'.');
% plot(lmarks(1,1 : n),lmarks(2,1 : n),'.','Markersize',20);

%% center the data around the mean of the landmarks
c = mean(lmarks,2);
lmarks = lmarks - c*ones(1,n);
k = 2; 
%% Step 2:  classical MDS on landmarks
B = lmarks'*lmarks;
[V,E] = eig(B);
[lam,isort] = sort(diag(E),'descend');
V = V(:,isort);
ind = find(lam > 0);
if isempty(ind)
    fprintf('The landmark distance squared matrix has no positive eigenvalues\n');
    return;
end
kplus = min(k,length(ind));
lamk = lam(1 : kplus);
L = diag(sqrt(lamk))*(V(:,1 : kplus)');
fprintf('Done with MDS on landmarks\n');
%% Step 3: distance-based triangulation
data = data - c*ones(1,Ndata);
dots = [lmarks,data];
% D = matrix of distances squared
fprintf('Computing D\n');
for i = 1 : n %N
    D(i,:) = sum((dots(:,i)*ones(1,N) - dots).^2,1);
end
fprintf('Done with computing D \n');
murow = zeros(n,1);
for i = 1 : n
    murow(i) = mean(D(i,:));
end
fprintf('Done with computing murow\n');
Lpseudo = diag((1./sqrt(lamk)))*(V(:,1 : kplus)');
X = zeros(n,N);
% size(Lpseudo)
% size(murow)
% size(D)
X = -0.5*Lpseudo*(D(1 : n,:) - murow*ones(1,N)); % matrix-vector multiplication
fprintf('Done with distance-based triangulation\n');
%% Step 4: PCA normalisation (optional). 
Xmean = mean(X,2);
Xhat = X - Xmean*ones(1,N);
[U,E] = eig(Xhat*Xhat');
[kap,isort] = sort(diag(E),'descend');
U = U(:,isort);
Xpca = U'*Xhat;
fprintf('Done with PCA normalization\n');
% figure(1); clf; hold on;
% plot(Xpca(1,:),Xpca(2,:),'.');
% plot(Xpca(1,1 : n),Xpca(2,1 : n),'.','Markersize',20);

newlmarks = Xpca(:,1 : n);
newdata = Xpca(:,n + 1 : N);
% figure(4); clf; hold on;
% plot(newdata(1,:),newdata(2,:),'.');
% plot(newlmarks(1,1 : n),newlmarks(2,1 : n),'.','Markersize',20);

% Output: k+ ? L? X? Xpca?
% actual number of output dimensions
% k+ -dimensional embedding of landmarks
% k+ -dimensional embedding of data
% PCA-normalised k+-dimensional embedding of data
end
