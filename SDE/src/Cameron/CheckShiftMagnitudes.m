function CheckShiftMagnitudes()
data = load('delta_net_CM2D.mat');
net = data.net; % coordinates of vertices of the delta-net
E = data.E; % edges of the delta-net
xgrid = data.xg;
ygrid = data.yg;
Ugrid = data.U;
dd = load('ab.mat’); % absolute value of the field b
ab = dd.ab;
mi = min(ab);
ma = max(ab);
figure;
axis square;
hold on; grid;
n = length(ab);
for i = 1 : n 
    c(i) =(ab(i) - mi)/(ma - mi); 
    plot3(net(1,i),net(2,i),ab(i),'.','Markersize',20,'color',[c(i),0,1-c(i)]);
    drawnow
end
plot3(0,0,0,'.','Markersize',40,'color',[c(i),0,1-c(i)]);
plot3(0,1,0,'.','Markersize',40,'color',[c(i),0,1-c(i)]);
plot3(1,1,0,'.','Markersize',40,'color',[c(i),0,1-c(i)]);



%%%%%%%
dd = load('LearnedSimulator_CM2D.mat');
b = dd.b;
ab = sqrt(sum(b.^2,1)); % compute absolute values of b
mi = min(ab);
ma = max(ab);
figure;
axis square;
hold on; grid;
n = length(ab);
for i = 1 : n 
    c(i) =(ab(i) - mi)/(ma - mi); 
    plot3(net(1,i),net(2,i),ab(i),'.','Markersize',20,'color',[c(i),0,1-c(i)]);
    drawnow
end
plot3(0,0,0,'.','Markersize',40,'color',[c(i),0,1-c(i)]);
plot3(0,1,0,'.','Markersize',40,'color',[c(i),0,1-c(i)]);
plot3(1,1,0,'.','Markersize',40,'color',[c(i),0,1-c(i)]);


end
