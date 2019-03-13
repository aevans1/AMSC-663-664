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


