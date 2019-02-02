function [x,j] = learned_simulator_step(x,i,new_S,neighbors,d,dt,delta)
%TODO: fix new_i, j stuff

%Given the learned simulator new_S, computes timestep at point x in chart i
%inputs:    x- initial point for simulator
%           i- chart index that x resides in
%           new_S- learned simulator, struct with values:
%                   T- transition mappings, a struct array
%                   with T(i,j).T represent the transition matrix from
%                   chart i to chart j
%                   B- diffusion coefficients for each chart, array of
%                   form B(:,n) for chart n
%                   C- local coordinates for net points of each chart,
%                   struct of form C(n,i).C: i-th neighbor of chart point
%                   y_n, represented in coordinates of chart n
%                   Sigma- drift coefficients for each chart, array of
%                   form C(:,n) for chart n
%                   mu- local coordinates for mean of common landmarks,
%                   structy array with mu(n,j).mu representing the mean
%                   landmark for y_j U y_n, represented in chart y_n
%                   Phi- 1-d case only, this is the MDS transformation
%                   multiplier, Phi(n).Phi = +-1, transformation for chart
%                   y_n
%           neighbors -struct array, each net point index has a struct with key value
%		             'nbr' that stores the indices of close by net points
%           d - intrinsic dimension of dynamics, problem specifics
%           delta - homogenization scale, affects density of sample net
%           dt = desired time-step length
%outputs:   x: new x value after timestep
%           j: new chart index for x

%read in Struct new_S
C = new_S.C;
B = new_S.B;
Sigma = new_S.Sigma;
T = new_S.T;
mu = new_S.mu;

net_nbr = neighbors(i).nbr;
num_nbr = length(net_nbr);
distances = [];
for n = 1:num_nbr
    j = net_nbr(n);
    distances(n) = norm(x - C(n,j).C);
end

%find global index for closest chart
new_i =net_nbr(find(distances == min(distances),1));

j = i; %j is new chart to be assigned to next sim step
if new_i ~= i
    
    %TESTING
    %         fprintf("i = %d     ",i);
    %         fprintf("neighbors = ");
    %         net_nbr
    %         fprintf("new_i = %d     \n",new_i);
    %         fprintf("T(%d,%d).T = ",i,new_i);
    %         T(i,new_i).T
    %         fprintf("mu(%d,%d).mu = ",i,new_i);
    %         mu(new_i,i).mu
    %
    %         fprintf("old x:");
    %         x
    %
    x = (T(i,new_i).T).' * (x - mu(i,new_i).mu) + mu(new_i,i).mu;
    
    %set new chart to i
    j = new_i;
end



%Forward Euler step
%take eta from gaussian N(0,eye(d))
eta = (randn(d));
x = x + B(:,j)*dt + eta*Sigma(:,j)*sqrt(dt);

%prevent escape from local chart
if norm(x) > (3/2)*delta
     x = (1/norm(x))*x*(2*delta - (1/2)*delta*exp(3 - (2/delta)*norm(x)));
end

end
