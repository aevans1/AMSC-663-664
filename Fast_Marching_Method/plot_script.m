x = linspace(0,1,257);
y = linspace(0,1,257);
[X,Y]=meshgrid(x,y);
U = load('U.txt');
V = load('V.txt');
contourf(X,Y,U - V);
