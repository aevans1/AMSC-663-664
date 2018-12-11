x = linspace(0,1,2049);
y = linspace(0,1,2049);
[X,Y]=meshgrid(x,y);
U = load('U.txt');
err = load('err.txt');
contourf(X,Y,U - err);

