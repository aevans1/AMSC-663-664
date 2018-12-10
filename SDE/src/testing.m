X = rand(5,10);
rho = @(x,y) norm(x-y);
d = 5;

L = X(:,1:3);

[new_L,eigvecs] = PCA(X,2)

test_X = (eigvecs.')*X;

[new_L,new_X] = LMDS(L,X,rho,2);
new_X - test_X




