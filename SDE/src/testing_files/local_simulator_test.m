%load relevant files before! list here
%TODO: ^^
%b = 3;
%Sigma = 2;
%rng(sum(clock));	
%p = 10000; %num_paths
%T = 1;
%dt = 0.01;
%N = floor(T/dt);
%Y = zeros(1,N+1); 
%
%Yzero = rand();
%path_end = zeros(1,p);
%for k = 1:p
%	Y(:,1) = Yzero;
%	Ytemp = Yzero;
%	%%%Compute SDE until time dt*N
%	for j = 1:N
%		eta = randn();
%		Ytemp = Ytemp + b*dt + Sigma*eta*sqrt(dt);
% 		%Y(:,j+1) = Y(:,j) + b*dt + sigma*eta*sqrt(dt);
%		Y(:,j+1) = Ytemp;
%	end
%	path_end(:,k) = Y(:,N+1);
%end
%	
%Y_mean = mean(path_end);
%true_mean = b*T + Yzero;
%Y_var = cov(path_end);
%true_var = Sigma^2*T;

B = new_S.B;
Sigma = new_S.Sigma;
rng(sum(clock));	
p = 10000; %num_paths
T = .2;
dt = t_0/5;
for n = 1:size(net,2)
	%Yzero = 1.6*rand() - 0.3;
	Yzero = 0;	
	N = floor(T/dt);
	Y = zeros(d,N+1); 
	total_Y = zeros(d,m*(N+1)); %collection of all m paths
  	for k = 1:p
		Y(:,1) = Yzero;
		Ytemp = Yzero;
		%%%Compute SDE until time dt*N
		for j = 1:N
			eta = randn(d,1);
			Ytemp = Ytemp + B(:,n)*dt + Sigma(:,:,n)*eta*sqrt(dt);
	 		Y(:,j+1) = Ytemp; 
		end
		total_Y(:,1 + (k-1)*(N+1): k*(N+1)) = Y;
		path_end(:,i) = Y(:,N+1);
    end

	Y_mean = mean(path_end);
	Y_var = cov(path_end);

	mean_error = abs((Y_mean - t_0*B(:,n))/(t_0*B(:,n)))
	var_error = abs((Y_var - Sigma(:,n)*Sigma(:,n).'*t_0)/(Sigma(:,n)*Sigma(:,n).'*t_0))
	mean_error = abs((Y_mean - B(:,n))/(B(:,n)))
	var_error = abs((Y_var - Sigma(:,n)*Sigma(:,n).')/(Sigma(:,n)*Sigma(:,n).'))

end




