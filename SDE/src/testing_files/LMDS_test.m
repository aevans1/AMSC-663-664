rng(1);

d = 1;
D = 1;

X = rand(D,5);
L = rand(D,2);
rho = @(x,y) norm(x - y);
Xshift = X - mean(L,2);
Lshift = L - mean(L,2);

[local_L,local_X] = LMDS(L,X,rho,d)
%local_X(:,1)./X(:,1)
%local_X(:,1)./Xshift(:,1)
%%local_L(:,1)./L(:,1)
%local_L(:,1)./Lshift(:,1)
%fprintf("next!\n");

[local_L,local_X] = LMDS(L,Xshift,rho,d);
%local_X(:,1)./X(:,1)
%local_X(:,1)./Xshift(:,1)
%local_L(:,1)./L(:,1)
%local_L(:,1)./Lshift(:,1)
%fprintf("\n");
%
[local_Lshift,local_Xshift] = LMDS(Lshift,Xshift,rho,d)
%local_Xshift(:,1)./X(:,1)
%local_Xshift(:,1)./Xshift(:,1)
%local_Lshift(:,1)./L(:,1)
%local_Lshift(:,1)./Lshift(:,1)
%fprintf("\n");
%[local_Lshift,local_Xshift] = LMDS(Lshift,X,rho,d);
%local_Xshift(:,1)./X(:,1)
%local_Xshift(:,1)./Xshift(:,1)
%local_Lshift(:,1)./L(:,1)
%local_Lshift(:,1)./Lshift(:,1)





