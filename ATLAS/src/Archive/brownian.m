%Computes Brownian motion over interva [0,T] with N steps
function W = brownian(T,N)
rng(sum(clock));
T = 1; N = 500; dt = T/N;
dW = sqrt(dt)*randn(1,N);   %brownian increments dW1,dW2,...
W = cumsum(dW);
end