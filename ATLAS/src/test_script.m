X = rand(5,1000);
L = X(:,1:500);
rho = @(x,y) norm(x - y);
LMDS(L,X,rho,2,false,true);