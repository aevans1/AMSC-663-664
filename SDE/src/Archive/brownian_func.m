T = 1;
N = 500;
dt = T/N;
t = [dt:dt:1]; %for plotting

M = 1000; %repeating M random walks
dW = sqrt(dt)*randn(M,N);
W = cumsum(dW,2); %each row is a random walk
U = my_func(W);
Umean = mean(U);
plot([0 t], [1, Umean] ); hold on;
plot([0 t], [ones(5,1) U(1:5,:)]);

function U = my_func(X)
    U = X.^2;
end
