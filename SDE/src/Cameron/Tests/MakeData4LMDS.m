function data = MakeData4LMDS(nu)
% nu = noise level
% nu = 0.0; 
% place points on a grid
nx = 30;
ny = 20;
x = linspace(-1.5,1.5,nx);
y = linspace(-1,1,ny);
[X,Y] = meshgrid(x,y);
dots = [X(:),Y(:)];
N = size(dots,1);
% compute the matrix of squared distances
D = zeros(N);
for j = 1 : N
    D(j,:) = sum((dots - ones(N,1)*dots(j,:)).^2,2);
end
% add noise
s = log(1 + nu); % nu = exp(s) - 1
noise = sqrt(s)*randn(N,N);
noise = 0.5*(noise + noise'); %symmetrize the noise;
data = D + noise;
end