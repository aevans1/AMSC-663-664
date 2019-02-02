%TODO: comments
close all;

x_init = [-1:0.01:2.5];
y_init = [-1:0.01:2];
[X,Y] = meshgrid(x_init,y_init);

figure;
U = arrayfun(@(X,Y) example_3(X,Y),X,Y);
W = (U < 4).*U;
[M,c] = contour(X,Y,W);
c.LineWidth = 2;
colormap(jet);
colorbar;
axis square;
hold on;

%2D example 

%remove mesh points with U <=10
reduced_X = (U < 10).*X;
reduced_Y = (U < 10).*Y;
init = [reduced_X(:)' ; reduced_Y(:)'];

rho = @(p_1,p_2) norm(p_1 - p_2);
delta = 0.2;

size(init)
[net,neighbors] = delta_net(init,delta,rho);
size(net)
for n = 1:size(net,2)
	nbrs = neighbors(n).nbr;
	num_nbr = size(nbrs,2);
	for i = 1:num_nbr
		nbr = nbrs(i);
		plot( [net(1,n),net(1,nbr)],[net(2,n),net(2,nbr)],'-o','Color','k');
		hold on;
	end
end

%hold on;

%%%%%%%%%%%%%%%%%%%
%TESTING: is the the correct function being plotted?
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


