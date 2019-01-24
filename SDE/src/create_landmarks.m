function landmarks = create_landmarks(S,delta_net,m,t_0)
%for each y in delta_net, run simulator m times for time t_0 each, keep
%path endpoints as landmarks corresponding to y
%inputs: S - SDE simulator
%		 delta-net - delta-net generated from delta_net function
%		 m - desired number of landmarks to generate for each net point
%		 t_0 - desired time of simulation for each call to S
%outputs: landmarks - array of land-marks, dim D x (m+1) x N
%			1st. dim is coordinate %of landmark point
%			2nd. dim is index of landmark point for given delta-net point x_n,
%		    	index 1 is the delta-net point x_n, next m are landmarks
%			3rd. dim is index of delta_net point
%e.g:
%A(:,:,n) is set of landmarks for delta-net point x_n
%A(:,1,:) is set of delta-net points
	
	D = size(delta_net,1); %assuming each column of delta_net is a data_point
	N = size(delta_net,2); %number of points in delta_net
	landmarks = zeros(D,m+1,N);
	for n = 1:N
	%generate m paths, save endpoints for delta_net point n
	%NOTE: %first point in landmark list is the %delta_net point itself
			landmarks(:,1,n) = delta_net(n);
			landmarks(:,2:m+1,n) = S(delta_net(n),m,t_0);
		end
end