function LMDStest()
nu = 0.1; % noise level
data = MakeData4LMDS(nu);
% data is a symmetricmatrix of distances squared
N = size(data,1); % data is N-by-N
k = 2; % the desired number of output dimensions;
n = 10; % the number of landmark points
%%
% Step 1: select landmark points (randomly)
P(1:N) = randperm(N)'; % random permutation 
lmarks = P(1:n); % landmarks
% Step 2:  classical MDS on landmarks
D = data(lmarks,lmarks); % distance squared matrix for landmarks
murow = mean(D,2); % row means of D
mu = mean(murow);
B = zeros(n);
for i = 1 : n
    for j = 1 : n
        B(i,j) = -0.5*(D(i,j) - murow(j) - murow(i) + mu);
    end
end
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
% Step 3: distance-based triangulation
Lpseudo = diag((1./sqrt(lamk)))*(V(:,1 : kplus)');
X = zeros(n,N);
X = -0.5*Lpseudo*(data(lmarks,:) - murow*ones(1,N)); % matrix-vector multiplication
% Step 4: PCA normalisation (optional). 
Xmean = mean(X,2);
Xhat = X - Xmean*ones(1,N);
[U,E] = eig(Xhat*Xhat');
[kap,isort] = sort(diag(E),'descend');
U = U(:,isort);
Xpca = U'*Xhat;

figure(1); clf; hold on;
plot(Xpca(1,:),Xpca(2,:),'.');
plot(Xpca(1,lmarks),Xpca(2,lmarks),'.','Markersize',20);

% Output: k+ ? L? X? Xpca?
% actual number of output dimensions
% k+ -dimensional embedding of landmarks
% k+ -dimensional embedding of data
% PCA-normalised k+-dimensional embedding of data
end
