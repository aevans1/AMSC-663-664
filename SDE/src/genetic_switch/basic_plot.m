figure(1); clf; hold on;
%contour(xgrid,ygrid,Ugrid,linspace(0,3.7,20),'Linewidth',1);
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
set(gca,'Fontsize',20);
for i = 1 : length(edges)
    xedge = [net(1,edges(i,1)),net(1,edges(i,2))];
	yedge = [net(2,edges(i,1)),net(2,edges(i,2))];
	zedge = [net(3,edges(i,1)),net(3,edges(i,2))];
    plot3(xedge,yedge,zedge,'color','k','LineWidth',0.25);
end
drawnow;


