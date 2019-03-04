close all;
load('atlas_2D');

x_init = [-1:0.01:2.5];
y_init = [-1:0.01:2];
[X,Y] = meshgrid(x_init,y_init);

figure;
U = arrayfun(@(X,Y) example_3(X,Y),X,Y);
W = (U < 4).*U;
[M,c] = contour(X,Y,W,20);
c.LineWidth = 1;
colormap(jet);
colorbar;
axis square;
hold on;
%%

%remove mesh points with U <=10
reduced_X = (U < 10).*X;
reduced_Y = (U < 10).*Y;
init = [reduced_X(:)' ; reduced_Y(:)'];
%%
%% Plot delta_net
for n = 1:size(net,2)
	nbrs = neighbors(n).nbr;
	num_nbr = size(nbrs,2);
	for i = 1:num_nbr
		nbr = nbrs(i);
		plot([net(1,n),net(1,nbr)],[net(2,n),net(2,nbr)],'-o','Color','k','MarkerFaceColor','b');
	end
end
hold off;

%
%
%%% Plot 1 Landmark per point
%for n = 1:size(net,2)
%	plotting_A = A(:,:,n);
%	for i = 2:3
%		plot( [net(1,n),plotting_A(1,i)],[net(2,n),plotting_A(2,i)],'--o','Color','r','MarkerFaceColor','r');
%		hold on;	
%	end
%	plot([net(1,n),net(1,n)],[net(2,n),net(2,n)],'-o','Color','k','MarkerFaceColor','b');
%
%end
%


%%%%%%%%%%%%%%%%%%%%%%
%Plot drift vectors for original v. ATLAS
%%%%%%%%%%%%%%%%%%%%%%
%%% plot drift vectors for original system
drift = [];
for n = 1:size(net,2)
	drift(:,n) = example_3_grad(net(1,n),net(2,n));
end
figure;
hold on; grid;


plot_num = size(net,2);
for n = 1:plot_num
	quiver(net(1,n),net(2,n),0.05*drift(1,n),0.05*drift(2,n),'AutoScale','off','MaxHeadSize',1.0,'color','k');
end
plot3(0,0,0,'.','Markersize',40,'color','b');
plot3(1.5,0,0,'.','Markersize',40,'color','b');
plot3(0.8,1.05,0,'.','Markersize',40,'color','b');

hold off;
%%% plot drift vectors for reduced ATLAS
figure;
hold on; grid;
for n = 1:plot_num
	quiver(net(1,n),net(2,n),0.05*B(1,n),0.05*B(2,n),'AutoScale','off','MaxHeadSize',1.0,'color','r');
end
plot3(0,0,0,'.','Markersize',40,'color','b');
plot3(1.5,0,0,'.','Markersize',40,'color','b');
plot3(0.8,1.05,0,'.','Markersize',40,'color','b');

hold off;

%%%plot both sets of drift vectors on same
figure;
hold on; grid;
plot_num = size(net,2);
for n = 1:plot_num
	quiver(net(1,n),net(2,n),0.05*drift(1,n),0.05*drift(2,n),'AutoScale','off','MaxHeadSize',1.0,'color','k');
	quiver(net(1,n),net(2,n),0.05*B(1,n),0.05*B(2,n),'AutoScale','off','MaxHeadSize',1.0,'color','r');
end
plot3(0,0,0,'.','Markersize',40,'color','b');
plot3(1.5,0,0,'.','Markersize',40,'color','b');
plot3(0.8,1.05,0,'.','Markersize',40,'color','b');

hold off;



%%%%%%%%%%%%%%%%%%%%%%
%Plot magnitudes of driver vectors for original v. ATLAS
%%%%%%%%%%%%%%%%%%%%%%
%%% find magnitudes of drift vectors for original system
mag_drift = [];
for n = 1:size(net,2)
	mag_drift(n) = norm(example_3_grad(net(1,n),net(2,n)));
end

min_drift = min(mag_drift);
max_drift = max(mag_drift);
figure;
hold on; grid;
for n = 1:size(net,2)
	color(n) = (mag_drift(n) - min_drift)/(max_drift - min_drift);
	plot3(net(1,n),net(2,n),mag_drift(n),'.','Markersize',20,'color',[color(n),0,1-color(n)]);
end
plot3(0,0,0,'.','Markersize',40,'color',[color(n),0,1-color(n)]);
plot3(1.5,0,0,'.','Markersize',40,'color',[color(n),0,1-color(n)]);
plot3(0.8,1.05,0,'.','Markersize',40,'color',[color(n),0,1-color(n)]);

hold off;
%% find magnitudes of drift vectors for ATLAS
mag_B = [];
for n = 1:size(net,2)
	mag_B(n) = norm(B(:,n));
end

min_B = min(mag_B);
max_B = max(mag_B);
figure;
hold on; grid;
for n = 1:size(net,2)
	color(n) = (mag_B(n) - min_B)/(max_B - min_B);
	plot3(net(1,n),net(2,n),mag_B(n),'.','Markersize',20,'color',[color(n),0,1-color(n)]);
end
plot3(0,0,0,'.','Markersize',40,'color',[color(n),0,1-color(n)]);
plot3(1.5,0,0,'.','Markersize',40,'color',[color(n),0,1-color(n)]);
plot3(0.8,1.05,0,'.','Markersize',40,'color',[color(n),0,1-color(n)]);
%
%%%%%%%%%%%%%%%%%%%
%TESTING: is the the correct function being plotted for the gradient?
%figure;
%example_3_norm_grad = @(a,b) (norm(example_3_grad(a,b))^2);
%Z = arrayfun(@(X,Y) example_3_norm_grad(X,Y),X,Y);
%W = (Z < 100).*Z;
%contour(X,Y,W);
%colorbar;
%axis square;
%%%%%%%%%%%%%%%%%%%
function out = example_3(x_1,x_2)
	x = [x_1 ; x_2];
	c = [1/5;1/5;1/6];
	p_1 = [0;0];
	p_2 = [1.5;0];
	p_3 = [0.8;1.05];
	out =  -log( exp( (1/c(1))*(-norm(x - p_1)^2)) + exp( (1/c(2))*(-norm(x - p_2)^2)) + exp( (1/c(3))*(-norm(x - p_3)^2)));
end
function out = example_3_grad(x_1,x_2)
	x = [x_1 ; x_2 ];
	c = [1/5;1/5;1/6];
	p_1 = [0;0];
	p_2 = [1.5;0];
	p_3 = [0.8;1.05];
	out =  -2*( ( (x - p_1)/c(1))*exp( (1/c(1))*(-norm(x - p_1)^2)) + ( (x - p_2)/c(2))*exp( (1/c(2))*(-norm(x - p_2)^2)) + ((x-p_3)/c(3))*exp( (1/c(3))*(-norm(x - p_3)^2)))/...
	( exp( (1/c(1))*(-norm(x - p_1)^2)) + exp( (1/c(2))*(-norm(x - p_2)^2)) + exp( (1/c(3))*(-norm(x - p_3)^2)));
end


