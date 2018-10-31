x = linspace(0,1,257);
y = linspace(0,1,257);
[X,Y]=meshgrid(x,y);
U = load('U.txt');
V = U;
%%
%load in different U from other git branch
U = load('U.txt');
%contourf(X,Y,U,20);
