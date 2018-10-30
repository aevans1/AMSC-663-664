x = linspace(0,1,1029);
y = linspace(0,1,1029);
[X,Y]=meshgrid(x,y);
U = load('U.txt');
contourf(X,Y,U,20);
