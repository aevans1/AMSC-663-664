function construction(params)
%Constructs delta-net, landmarks, ATLAS, and local SDE simulator

%inputs:params struct, which contains
%		S - SDE simulator
%       delta - homogenization scale, affects density of sample net
%		rho -  given distance function
%		m - desired number of landmarks to generate for each net point
%       p - num sample paths for each point in net
%		t0 - desired time of simulation for each call to S
%       d - intrinsic dimension of dynamics, problem specific
%		net_info - struct for delta_net, contains
%			(N = number of net points)
%			net - d X N array, columns are data points of delta net
%			neighbors - N x max_deg array, row j is neighbor indices of net	point j
%			deg - N x 1 vector, entry j is number of neighbors(degree) of net point j
%			max_deg - maximum enty of deg
%outputs: new_S - struct with key values that specify parameters for
%                   simulating a constant coefficient SDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read in data, incoporate flags and extra parameters

%%%Allow parameter inputs to construcion or load up saved inputs
if ~exist('params','var')
    try
        fprintf("No input parameters to construction, attempting to load saved params in current_driver.mat \n");
        data = load('current_driver.mat');
    catch
        fprintf("Error: No saved params under current_driver.mat, please generate data parameters in atlas_driver.m \n");
        return;
    end
    params = data.params;
end

%%%Check that param struct has non-empty fields
fields = {'S','delta','rho','m','p','t0','d','net_info'};

for i = 1 : length(fields)
    if ~isfield(params,fields{i})
        fprintf("Error: params.%s is not set, please re-initialize param struct with %s value \n", fields{i},fields{i});
        return;
    end
    if isempty(getfield(params,fields{i}))
        fprintf("Error: params.%s is empty, please re-initialize param struct with %s value \n", fields{i},fields{i});
        return;
    end
end

%%%Read in struct parameters
S = params.S;
delta = params.delta;
rho = params.rho;
m = params.m;
p = params.p;
t0 = params.t0;
d = params.d;

net_info= params.net_info;
net = net_info.net;
neighbors = net_info.neighbors;
edges = net_info.edges;
deg = net_info.deg;
max_deg = net_info.max_deg;

%%% Flags, extra parameters go here
timestamp_save = false; %if true, saves a mat file with timestamp
testing = true; %if true, runs various test blocks of code in file
vary_t0 = false; %if true, allow for varying t0 values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Initialize arrays
D = size(net,1); %each column of delta_net is a data_point, D is ambient dimension of space
N = size(net,2); %number of points in delta_net

X = zeros(D,p,N); 							%simulated points for learning ATLAS
A = zeros(D,m+1,N); 						%landmark set for each delta net point
L = zeros(D,(max_deg + 1)*(m+1),N); 		%landmarks + neighbors landmarks for each delta net point
embed_L = zeros(d,(max_deg + 1)*(m+1),N);   %mapping of L via LMDS
embed_X = zeros(d,p,N); 					%mapping of X via LMDS
Phi = zeros(d,N); 							%for d = D = 1 case, gives LMDS mapping coefficient
T = zeros(d,d,N,N); 						%T is the collection of transition maps from chart to chart
c = zeros(d,N,N); 							%c is the collection of chart centers for each chart
b = zeros(d,N); 							% collection of drift coeff. for local SDEs
sigma = zeros(d,d,N); 						% collection of diffusion matrices for local SDEs;
mu = zeros(d,N,N); 							%mu is the collection of average landmark for each chart


%%%TESTING: using experimentally verified t0_vals, only valid for current
%%%gswitch net
if vary_t0
   load('t0_list.mat','t0_list'); 
end

%%%Create landmarks
for n = 1:N
    %generate m paths, save endpoints for delta_net point n
    %NOTE: %first point in landmark list is the %delta_net point itself
    A(:,1,n) = net(:,n);
    
    %use default t0 unless we're in 'vary t0' mode, then use from list
    if ~vary_t0
        tval = t0;
    else
        tval = t0_list(n);
    end
    
    
    A(:,2:m+1,n) = S(net(:,n),m,tval);
end

%%%Create atlas
for n = 1:N
    fprintf("Constructing for net point %d \n",n);
    
    
    %%%Step 1:Create landmarks and simulator trajectories(X) for net point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%Simulate p paths around net point y_n
    if ~vary_t0
        tval = t0;
    else
        tval = t0_list(n);
    end
    X(:,1:p,n) = S(net(:,n),p,tval);

    %TODO: remove?
    %TESTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Optional: run tests per net point n, see comments in test functions
%     if testing
%         mean_net_dist = test_mean_net_dist(n,D,net,delta,t0,X);
%         mean_net_distances(n) = mean_net_dist;
%         
%         new_t0 = t0; %initialize new t0 value at original
%         
%         save('atlas_diagnostics.mat','mean_net_distances','t0_list','n','N');
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
    Xn = X(:,1:p,n);
    
    
    %%%%Take union of all neighboring landmarks to y_n
    %%%%store points from neighboring landmarks as collection of columns %associated to y_n
    %local_nbr contains global indices of neighbor net points to netpoint n
    num_nbr = deg(n);
    local_nbr = neighbors(n,1:num_nbr);
    L(:,1:(m+1)*(num_nbr + 1),n) = reshape(A(:,:, [n local_nbr]),D,(m+1)*(num_nbr + 1));
    Ln = L(:,1:(m+1)*(num_nbr + 1),n);
    
    %TESTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Optional: run tests per net point n, see comments in test functions
    if testing
        test_landmark_index(n,m,N,neighbors,deg,A,Ln);
    end
    
    %%%Step 2:Create chart for net point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [embed_Ln,embed_Xn] = LMDS(Ln,Xn,rho,d);
    
    %%%For d = D == 1 case, track what Phi is, the LMDS mapping
    local_Phi = 0;
    if d ==1 && D ==1
        %Here Phi is whatever was multiplied by X's to get 'embedding'
        temp_Xn = Xn - mean([Ln Xn],2);
        local_Phi = embed_Xn(:,1)/temp_Xn(:,1);
        %local_Phi = embed_Xn(:,1)/Xn(:,1);
        
    end
    Phi(d,n)= local_Phi;
    
    %%%translate all charts to be centered around embedded net point n
    chart_center = embed_Ln(:,1);
    embed_Ln = embed_Ln - chart_center;
    embed_Xn = embed_Xn - chart_center;
    
    embed_L(:,1:(m+1)*(num_nbr + 1),n) = embed_Ln;
    embed_X(:,:,n) = embed_Xn;
    
    
    %%%Step 3:Construct local constant coefficient SDE simulator at net point
    %%%		dXt = b_n dXt + sigma_n dW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [T,c,b,sigma,mu] = construct_local_SDE(n,net_info,embed_L,embed_X,p,tval,m,d,T,c,b,sigma,mu);
    
end
%End loop through all net points

%%%Save new simulator parameters for ATLAS
new_S.T = T;
new_S.b = b;
new_S.c = c;
new_S.sigma = sigma;
new_S.mu = mu;
new_S.Phi = Phi;


save('current_atlas.mat','new_S','L','embed_L');


%TODO:remove?
% 
% if testing
%     %%%Save data collected on simulation distances and new t0 values
%     atlas_diagnostics.mean_net_distances = mean_net_distances;
%     save('current_atlas.mat','atlas_diagnostics','-append');
%     
%     if vary_t0
%        save('current_atlas.mat','t0_list','-append');
%     end
%     
% end

if timestamp_save
    filename =[datestr(now, 'dd_mmm_yyyy_HH_MM'),'_','d',num2str(d),'_','construction'];
    save(filename);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,c,b,sigma,mu] = construct_local_SDE(n,net_info,embed_L,embed_X,p,t0,m,d,T,c,b,sigma,mu)
%inputs:
%		n - current net point, function constructs local_SDE for chart n
%		net_info - struct for delta_net, contains
%			(N = number of net points)
%			net - d X N array, columns are data points of delta net
%			neighbors - N x max_deg array, row j is neighbor indices of net	point j
%			deg - N x 1 vector, entry j is number of neighbors(degree) of net point j
%			max_deg - maximum enty of deg
%		embed_L - LMDS applied to L_n, d x p array
%		embed_X - LMDS applied to X_n, d x num_nbr array
%		p - num sample paths for each point in net
%		t0 - desired time of simulation for each call to S
%		m - desired number of landmarks to generate for each net point
%		d - intrinsic dimension of dynamics, problem specifics
%		T - d X d X N X N array containing transition mappings between charts,
%			updated each iteration of construct_local_SDE,
%		    T(:,:,i,j) is a d x d matrix which shifts chart coords from chart i to j
%		c - d X N X N array containing local coordinates of neighboring charts,
%			 updated each iteration of construct_local_SDE,
%			 c(:,i,j) is a d x 1 vector, isj-th neighbor of net point n,
%		     expressed in chart n coords(LMDS wrt n's landmarks)
%		b - d X 1 vector, drift coordinate for n's local SDE
%		sigma - d x d x N array, sigma(:,:,n) is d x d array, diffusion coords
%		for n's local SDE, updated each iteration of construct_local_SDE
%		mu - d X N X N array containing mean landmarks of local neighboring charts,
%		     updated each iteration of construct_local_SDE,
%		    mu(:,i,j) is a d x 1 vector, the average of chart j's landmarks
%		    expressed in chart n coords
%outputs:
%		T,c,mu,b,sigma,mu: same as above, updated for net point n

%Note: following notation in ATLAS paper, 'y_n' is the embedded net point n

%%%Read in net_info and chart t embeddings
%Note: embed_Ln = [y_n |y_n landmarks |nbr1 |nbr1 landmarks |nbr2|nbr2 landmarks...];

neighbors = net_info.neighbors;
deg = net_info.deg;
max_deg = net_info.max_deg;
num_nbr_self = deg(n) + 1;
net_nbr = [n neighbors(n,1:deg(n))];

embed_Ln = embed_L(:,1:(m+1)*(num_nbr_self),n);
embed_Xn = embed_X(:,:,n);

%%%Compute chart c_n for y_n and center around embedded y_n
%NOTE: indexing to find nbr is based on that
%net_nbr = [n nbr1 nbr2 ]
%embed_Ln = [y_n landmarks| y_nbr1 landmarks |y_nbr2 ...]
for i = 1: num_nbr_self
    c(:,n,net_nbr(i)) = embed_Ln(:,(i-1)*(m+1) + 1);
end

%%%Compute diffusion, drift coefficients around y_n
b(:,n) = (1/ (p*t0) )*sum(embed_Xn,2);%drift for y_n
sigma(:,:,n) = sqrt(1/t0)*sqrtm(cov((embed_Xn).')); %diffusion for y_n

%%%Compute switching maps
%NOTE: indexing to find nbr is based on that
%net_nbr = [n nbr1 nbr2 ]
%embed_Ln = [y_n landmarks| y_nbr1 landmarks |y_nbr2 ...]
for i = 1 : num_nbr_self
    if net_nbr(i) < n
        j = net_nbr(i);	%global index of neighbor of n
        
        %%%Find embeddings of chart overlaps between chart for net point n and nbr net point j
        %L_jn is Phi324_n(A_n U A_j), the embedding of landmarks for net points n,j in chart n's coordinates
        %L_jn is Phi_j(A_n U A_j), the embedding of landmarks for net points n,j in chart j's coordinates
        
        L_nj = zeros(d,2*m + 2);
        L_jn = zeros(d,2*m + 2);
        
        L_nj(:,1:m+1) = embed_Ln(:,1:m+1);  %A_n wrt y_n embedding
        L_nj(:,m+2:2*m+2) = embed_Ln(:, (i - 1)*(m+1) + 1: i*(m+1)); %A_j wrt y_n embedding
        mu_nj = mean(L_nj,2);
        
        embed_Lj = embed_L(:,1:(m+1)*(deg(j) + 1),j); %embeddings wrt nbr of y_n
        nbr_switch_ind = [j neighbors(j, 1 : deg(j))];
        ind = find(nbr_switch_ind ==n);
        
        L_jn(:,1:m+1) = embed_Lj(:,(ind-1)*(m+1) + 1: ind*(m+1)); %A_n  wrt nbr j embedding
        L_jn(:,m+2:2*m+2) = embed_Lj(:,1:m+1); %A_j wrt nbr j embedding
        mu_jn = mean(L_jn,2);
        
        %%%Compute transition mapping T_nj between the two charts
        %T_nj changes x from chart n coords to chart j coords
        T_nj = (L_jn - mu_jn)*pinv(L_nj - mu_nj);
        T_jn = (L_nj - mu_nj)*pinv(L_jn - mu_jn);
        
        mu(:,n,j) = mu_nj;
        mu(:,j,n) = mu_jn;
        T(:,:,n,j) = T_nj;
        T(:,:,j,n)= T_jn;
        %end transition map computing for neighbor
    end
    %end computation of transition maps
end
end

function test_landmark_index(n,m,N,neighbors,deg,A,Ln)
%Checking that Ln data struct is properly organized
num_nbr = deg(n);
local_nbr = neighbors(n,1:num_nbr);
%fprintf("TESTING: index check for landmarks \n");
for i = 1:num_nbr
    ind = local_nbr(i);
    if(Ln(:,i*(m+1) + 1:(i+1)*(m+1)) ~= A(:,:,ind))
        fprintf("index error at landmark collection \n");
    end
end
end

function [dist_to_net] = test_mean_net_dist(n,D,net,delta,t0,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO: replace euclidean norm with more general norm
%t0 simulation time should be chosen so that
%E[X - net(:,n)] ~ delta, likewise for landmarks L
if D == 1
    dist_to_net =  mean(abs(X(:,:,n) - net(:,n)),2);
else
    dist_to_net =  mean(vecnorm(X(:,:,n) - net(:,n)),2);
end
fprintf("avg distance from X to net point n:	%f \n",dist_to_net);
fprintf("	delta: 	%f\n",delta);
fprintf("	t_0: 	%f\n",t0);

% %%%Recompute distance to net, output values
% if D == 1
% 	dist_to_net =  mean(abs(X(:,:,n) - net(:,n)),2);
% else
% 	dist_to_net =  mean(vecnorm(X(:,:,n) - net(:,n)),2);
% end
% fprintf("avg distance from X to net point n:	%f",dist_to_net);
% fprintf("new t0: %f \n",new_t0);

end
