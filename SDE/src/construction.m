function [new_Sim,neighbors,net] = construction(S,init,delta,rho,m,p,t_0,d,load_net)
%Constructs delta-net, landmarks, ATLAS, and SDE simulator
%inputs: S - SDE simulator
%        init - D X N matrix, set of N vectors in R^d
%        delta - homogenization scale, affects density of sample net
%		 rho -  given distance function
%		 m - desired number of landmarks to generate for each net point
%        p - num sample paths for each point in net
%		 t_0 - desired time of simulation for each call to S
%        d - intrinsic dimension of dynamics, problem specifics
%        load_net - boolean, true value loads up atlas variables specified in
%        			'preload.m', such as a delta net and neighbors
%outputs: new_Sim - struct with key values that specify parameters for
%                   simulating a constant coefficient SDE
%         net - columns are data points, all distances >= delta apart
%  		  neighbors - struct array, each net point index has a struct with key value
%		 'nbr' that stores the indices of close by net points


%if true for d = D, no LMDS computed and drift diffusion are from the
%original problem
no_LMDS = true;


%%%Create delta_net from initial points, create landmarks for delta_net
%	TODO: better way to do this? Loading up structs if saved
if ~exist('load_net','var')
    load_net = false;
end

if load_net == false
    [net,neighbors] = delta_net(init,delta,rho);
else
    fprintf("Using loaded variables for delta net and neighbors\n");
    load('preload','net','neighbors');
    fprintf('net is...\n');
    net
end

%fprintf("delta net is.... \n");
%net

D = size(net,1); %assuming each column of delta_net is a data_point
N = size(net,2); %number of points in delta_net

%%%Initialize arrays, and create landmarks from the net
X = zeros(D,p,N);
L(N).L = {};
embed_L(N).L = {}; %LMDS embeddings of landmarks
embed_X(N).X = {}; %LMDS embeddings of simulating points
Phi(N).Phi =[]; 	%multipliers for 1-dim chart mappings

T(N,N).T = {}; %T is the collection of transition maps from chart to chart
C(N,N).C = {}; %C is the collection of chart centers for each chart
mu(N,N).mu = {}; %mu is the collection of average landmark for each chart

%%%Create landmarks
A = zeros(D,m+1,N); %landmark set
for n = 1:N
    %generate m paths, save endpoints for delta_net point n
    %NOTE: %first point in landmark list is the %delta_net point itself
    A(:,1,n) = net(:,n);
    A(:,2:m+1,n) = S(net(:,n),m,t_0);
end

%%%Create atlas
for n = 1:N
    fprintf("Constructing for net point %d \n",n);
    
    %simulate p paths around net point y_n
    X(:,1:p,n) = S(net(:,n),p,t_0);
    
    %%%%Take union of all neighboring landmarks to y_n
    %net_nbr contains global indices of neighbor net points to netpoint n
    net_nbr = neighbors(n).nbr;
    num_nbr = length(net_nbr);
    
    %store points from neighboring landmarks as collection of columns
    %associated to y_n
    L(n).L = reshape(A(:,:,net_nbr),D,(m+1)*num_nbr);
    
    %TESTING: index check
    fprintf("TESTING: index check for landmarks \n");
    temp = L(n).L;
    for i = 1:num_nbr
        ind = net_nbr(i);
        if(temp(:,(i-1)*(m+1) + 1:i*(m+1)) ~= A(:,:,ind))
            fprintf("index error at landmark collection \n");
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TESTING: distances to net points from trajectories
    %according to ATLAS paper, t_0 should be chosen so that
    %E[X - net(:,n)] ~ delta, likewise for landmarks L
    %if D == 1
    %	dist_to_net =  mean(abs(X(:,:,n) - net(:,n)),2);
    %else
    %	dist_to_net =  mean(vecnorm(X(:,:,n) - net(:,n)),2);
    %end
    
    %fprintf("TESTING: avg distance from X to net point n:	%f",dist_to_net);
    %fprintf("	delta: 	%f\n",delta);
    %[M,I] = max(abs(X(:,:,n) - net(:,n)));
    %fprintf("TESTING: max distance from X to net point n:	%f, index %d \n",M,I);
    
    %%%same test as above, but for landmarks
    %if D == 1
    %	dist_to_net =  mean(abs(L(n).L - net(:,n)),2);
    %else
    %	dist_to_net =  mean(vecnorm(L(n).L - net(:,n)),2);
    %end
    %fprintf("TESTING: avg distance from landmarks L to net point n:	%f \n",dist_to_net);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%Step 1:Create chart for each net point
    
    if no_LMDS
        fprintf("TEMPORARY: no LMDS\n");
        local_L = L(n).L;
        local_X = X(:,:,n);
        local_Phi = 1;
    else
        [local_L,local_X,local_Phi] = create_chart(X(:,:,n),L(n).L,rho,d,D);
    end
    
    embed_L(n).L = local_L;
    embed_X(n).X = local_X;
    Phi(n).Phi = local_Phi;
    
    %TESTING: keeping track of original centers of each chart, i.e
    %images of net points, THEN centering
    embed_L(n).center = local_L(:,1);
    
	
	%TESTING: for NO_LMDS, don't shift by means
	if no_LMDS
		embed_X(n).X = local_X;
	else
		embed_L(n).L = embed_L(n).L - embed_L(n).center;
    	embed_X(n).X = local_X - embed_L(n).center;
   	end 
    %TESTING: redundant, but tracking these in the simulator struct as
    %well
    centers(:,n) = embed_L(n).center;
    
    %%%Step 2: constructing local SDE
    [T,C,B(:,n),Sigma(:,:,n),mu] = construct_local_SDE(n,neighbors,embed_L,embed_X,p,t_0,m,d,T,C,mu,no_LMDS);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TESTING: distances to net points from trajectories
    %if D ==1
    %	dist_to_net =  mean(abs(local_X),2);
    %else
    %	dist_to_net =  mean(vecnorm(local_X),2);
    %end
    %fprintf("TESTING: avg distance from projected X to projected net point n:	%f",dist_to_net);
    %fprintf("	2delta: 	%f\n",2*delta);
    %[M,I] = max(abs(local_X));
    %fprintf("TESTING: max distance from projected X to projected net point n:	%f, index %d \n",M,I);
    
    %if d ==1
    %	dist_to_net =  mean(abs(local_L),2);
    %else
    %	dist_to_net =  mean(vecnorm(local_L),2);
    %end
    %fprintf("TESTING: avg distance from projected L to projected net point n:	%f",dist_to_net);
    %fprintf("	2delta: 	%f\n",2*delta);
    %[M,I] = max(abs(local_L));
    %fprintf("TESTING: max distance from projected L to projected net point n:	%f, index %d\n",M,I);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
new_Sim.T = T;
new_Sim.B = B;
new_Sim.C = C;
new_Sim.Sigma = Sigma;
new_Sim.mu = mu;

%TESTING: tracking original chart centers for use in simulation
new_Sim.centers = centers;


%matters for D = d = 1 case only!
new_Sim.Phi = Phi;

filename =[datestr(now, 'dd_mmm_yyyy_HH_MM'),'_','d',num2str(d),'_','construction'];

save(filename);
save('current_construction');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [local_L,local_X,Phi] = create_chart(X,L,rho,d,D)
%TODO:comments

%tic
[local_L,local_X] = LMDS(L,X,rho,d);
%LMDS_time = toc

Phi = [];
if d ==1 && D ==1
    %%%%Note:D = d = 1 case only!
    %Here Phi is whatever was multiplied by X's to get embedding
    temp_X = X - mean(L,2);
    Phi = local_X(:,1)/temp_X(:,1);
end
%%%%%

%TODO: likely remove this
%%%Centering all data for chart around chart center
% 		center = local_L(:,1);
% 		local_L = local_L - center;
% 		local_X = local_X - center;

end

function [T,C,B,Sigma,mu] = construct_local_SDE(n,neighbors,embed_L,embed_X,p,t_0,m,d,T,C,mu,no_LMDS)
%TODO: comments
if ~exist('no_LMDS') no_LMDS = false; end


%net_nbr contains global indices of neighbor net points to netpoint n
net_nbr = neighbors(n).nbr;
num_nbr = length(net_nbr);

local_L = embed_L(n).L;
local_X = embed_X(n).X;

%select neighboring delta-net points:
%local_L = [y_n |y_n landmarks |nbr1 |nbr 1landmarks |nbr2|nbr2
%landmarks...];
ind = [1:m+1:(m+1)*num_nbr];

%%%Compute chart 'C_n' for y_n and center around embedded y_n
for i = 1:num_nbr
    j = net_nbr(i); %global index of neighbor
    C(n,j).C = local_L(:,ind(i)); %find nbrs locations relative to chart point n
end

%%%Compute Diffusion coefficients, drift coefficients around y_n
if no_LMDS
    fprintf("TEMPORARY: drift obtained from potential, diffusion is Identity \n");
    %TESTING: setting drift from original problem

    %TODO: line here selecting based off of example
    %B = example_1_grad(embed_L(n).center);
    B = example_3_grad(embed_L(n).center);
    
    %TESTING: setting diffusion from identity
    %option for true diffusion #1:
    Sigma = eye(d);
else
    B = (1/ (p*t_0) )*sum(local_X,2); %drift for y_n
    Sigma = sqrt(1/t_0)*sqrtm(cov(local_X.')); %diffusion for y_n
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%TESTING: index checks
%     fprintf("neighbors of net(:,%d): \n", n);
%     net_nbr

%%%Compute switching maps
%NOTE: net_nbr = [n idx1 idx 2 ...], idx1 < idx 2 < ...
for i = 1:num_nbr
    if net_nbr(i) < n
        j = net_nbr(i);	%global index of neighbor of n
        %TESTING:index checks
        %fprintf("neighbor of net(:,%d), global index: %d, local_index: %d \n",n,j,i);
                
        %L_j' in atlas paper
        local_nbr_L = embed_L(j).L; %embeddings wrt nbr of y_n
        
        %NOTE: indexing to find nbr is based on that 
        %net_nbr = [n nbr1 nbr2 ] 
        %[y_n landmarks| y_nbr1 landmarks ...]
        
        %NOTE: From paper: this is L_k,j
        %embedded landmarks pair wrt y_n
        L_nj(:,1:m+1) = local_L(:,1:m+1);  %y_k wrt y_k embedding
        L_nj(:,m+2:2*m+2) = local_L(:, (i - 1)*(m+1) + 1: i*(m+1)); %nbr wrt y_k embedding
        mu(n,j).mu = mean(L_nj,2);
        
        %%%find y_n wrt nbr embedding
        nbr_switch_ind = neighbors(j).nbr;
        ind = find(nbr_switch_ind ==n);
        
        %From paper: this is L_j,k
        %L_jn(:,1:m+1) = local_nbr_L(:,1:m+1);  %nbr i wrt nbr i embedding
        %L_jn(:,m+2:2*m+2) = local_nbr_L(:,(ind-1)*(m+1) + 1: ind*(m+1)); %y_n  wrt nbr i embedding
        
        %TESTING: swapping around the L_jn to match L_nj
        L_jn(:,1:m+1) = local_nbr_L(:,(ind-1)*(m+1) + 1: ind*(m+1)); %y_n  wrt nbr i embedding
        L_jn(:,m+2:2*m+2) = local_nbr_L(:,1:m+1); %nbr i wrt nbr i embedding
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TESTING: for test case with no LMDS, L_jn should equal L_nj
        %if ~all(sort(L_nj +embed_L(n).center) == sort(L_jn + embed_L(j).center))
        if ~all(L_nj +embed_L(n).center == L_jn + embed_L(j).center)
            fprintf("landmark pairs not equal!\n");
            L_nj + embed_L(n).center
            L_jn + embed_L(j).center
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        mu(j,n).mu = mean(L_jn,2);
        
        %T(n,j).T is transition map from chart n to chart j
        %i.e, T(n,j).T*x changes x from chart n coords to chart j
        %coords
        T(n,j).T = (L_jn - mu(j,n).mu)*pinv(L_nj - mu(n,j).mu);
        T(j,n).T = (L_nj - mu(n,j).mu)*pinv(L_jn - mu(j,n).mu);
        
    %end transition map computing for neighbor
    end
%end computation of transition maps

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
%Example 1: 1d smooth potential
%%%%%%%%
function out = example_1(x)
%example, 1d smooth potential
out = 16*x^2*(x-1)^2;
end
function out = example_1_grad(x)
%%% -grad(example_1(x))
out = -32*x*(x-1)*(2*x-1);
end

%%%%%%%%
%Example 2: 1D rough Potential
%%%%%%%
function out = example_2(x)
out = 16*x^2*(x-1)^2 + (1/6)*cos(100*pi*x);
end

function out = example_2_grad(x)
% -grad(example_2(x))
out = -32*x*(x-1)*(2*x-1) + (1/6)*(100*pi)*sin(100*pi*x);
end

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
