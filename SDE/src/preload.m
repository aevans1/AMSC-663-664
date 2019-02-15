%TODO: comments
%put any code here for computing delta nets or other parameters needed to use
%on multiple atlas runs


%%%ATLAS Example 3 delta net code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 3: Smooth 2-D Potential from ATLAS paper, ex. 5.3.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = @(x,y) norm(x - y);
delta = 0.2; 

%initial  point set for generating delta-net, [-1,2.5] x [-1,2] grid
%each col of init is (x,y) vector
init_x = [-1:0.01:2.5];
init_y = [-1:0.01:2]; 
[init_X,init_Y] = meshgrid(init_x,init_y);
temp_init = [init_X(:)';init_Y(:)'];
Z = [];
for i = 1:size(temp_init,2)
	Z(i) = example_3(temp_init(:,i));
end

%Removing all grid points with potential value below 10;
threshold = [Z < 10; Z < 10];
init = threshold.*temp_init;
init = init(:,any(init,1)); %removing zero columns

[net,neighbors] = delta_net(init,delta,rho);

save('preload');
%%%%%%%%
%Example 3: 2d Smooth Potential
%%%%%%%%
function out = example_3(x)
	%input: x is a vector in R^2
	%output: out is a scalar
	c = [1/5;1/5;1/6];
	p_1 = [0;0];
	p_2 = [1.5;0];
	p_3 = [0.8;1.05];
	out =  -log( exp( (1/c(1))*(-norm(x - p_1)^2)) + ...
		exp( (1/c(2))*(-norm(x - p_2)^2)) + ...
		exp( (1/c(3))*(-norm(x - p_3)^2)));
end

function out = example_3_grad(x)
	%-grad(example_3(x))
	%input: x is a vector in R^2
	%output: out is a vector in R^2
	
	c = [1/5;1/5;1/6];
	p_1 = [0;0];
	p_2 = [1.5;0];
	p_3 = [0.8;1.05];
	out =  -2*( ( (x - p_1)/c(1))*exp( (1/c(1))*(-norm(x - p_1)^2)) + ...
		( (x - p_2)/c(2))*exp( (1/c(2))*(-norm(x - p_2)^2)) +  ...
		((x-p_3)/c(3))*exp( (1/c(3))*(-norm(x - p_3)^2)))/...
		(exp( (1/c(1))*(-norm(x - p_1)^2)) + ...
		exp( (1/c(2))*(-norm(x - p_2)^2)) + ...
		exp( (1/c(3))*(-norm(x - p_3)^2)));
end



