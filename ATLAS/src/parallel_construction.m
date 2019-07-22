function parallel_construction(params)
%NOTE: identical to 'construction.m', but with some parellelization

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
%% Read in data, incoporate flags and extra parameters

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
testing = false; %if true, runs various test blocks of code in file
vary_t0 = false; %if true, allow for varying t0 values

% switches metric to 2-norm, allows for faster computations
if ~isfield(params,'euclidean')
    euclidean = false;
else
    euclidean = params.euclidean;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main Code

%%%Initialize arrays
D = size(net,1); %each column of delta_net is a data_point, D is ambient dimension of space
N = size(net,2); %number of points in delta_net
A = zeros(D,m+1,N); 						%landmark set for each delta net point
L = zeros(D,(max_deg + 1)*(m+1),N); 		%landmarks + neighbors landmarks for each delta net point
embed_L = zeros(d,(max_deg + 1)*(m+1),N);   %mapping of L via LMDS
embed_X = zeros(d,p,N); 					%mapping of X via LMDS
Phi = zeros(d,N); 							%for d = D = 1 case, gives LMDS mapping coefficient

%%%TESTING%%%%%%%%%%%%%%%%%%%
%Option to use varying t0_val chart by chart, current mat file only valid for current
%%%     gswitch net
if vary_t0
    load('t0_list.mat','t0_list');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Create landmarks
parfor n = 1:N
    
    %use default t0 unless we're in 'vary t0' mode, then use from list
    if ~vary_t0
        tval = t0;
    else
        tval = t0_list(n);
    end
    
    %generate m paths, save endpoints for delta_net point n
    %NOTE: %first point in landmark list is the %delta_net point itself
    %    A(:,1,n) = net(:,n);
    A(:,:,n) = [net(:,n) S(net(:,n),m,tval)];
end

%%%Create atlas
parfor n = 1:N
    fprintf("Constructing for net point %d \n",n);
    
    %%%Step 1:Create landmarks and simulator trajectories(X) for net point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Simulate p paths around net point y_n
    if ~vary_t0
        tval = t0;
    else
        tval = t0_list(n);
    end
    Xn = S(net(:,n),p,tval);
        
    %TESTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Optional: run tests per net point n, see comments in test functions
    if testing
        mean_net_dist = diagn_mean_net_dist(n,D,net,delta,t0,X);
        mean_net_distances(n) = mean_net_dist;
        
        new_t0 = t0; %initialize new t0 value at original
        
        save('atlas_diagnostics.mat','mean_net_distances','t0_list','n','N');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%Take union of all neighboring landmarks to y_n
    %   store points from neighboring landmarks as collection of columns %associated to y_n
    %   local_nbr contains global indices of neighbor net points to netpoint n
    num_nbr = deg(n);
    local_nbr = neighbors(n,1:num_nbr);
    Ln = reshape(A(:,:, [n local_nbr]),D,(m+1)*(num_nbr + 1));
    L(:,:,n) = [Ln zeros(D, (max_deg + 1)*(m+1) - (m+1)*(num_nbr + 1))];
    
    %TESTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%Run tests per net point n, see comments in test functions
    if testing
        test_landmark_index(n,m,N,neighbors,deg,A,Ln);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%Step 2:Create chart for net point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [embed_Ln,embed_Xn] = LMDS(Ln,Xn,rho,d,true,euclidean);
    
    %%%For d = D == 1 case, track what Phi is, the LMDS mapping
    local_Phi = 0;
    if d ==1 && D ==1
        %Here Phi is whatever was multiplied by X's to get 'embedding'
        temp_Xn = Xn - mean([Ln Xn],2);
        local_Phi = embed_Xn(:,1)/temp_Xn(:,1);
    end
    Phi(d,n)= local_Phi;
    
    %%%translate all charts to be centered around embedded net point n
    chart_center = embed_Ln(:,1);
    embed_Ln = embed_Ln - chart_center;
    embed_Xn = embed_Xn - chart_center;
    
    embed_L(:,:,n) = [embed_Ln zeros(d, (max_deg + 1)*(m+1) - (m+1)*(num_nbr + 1))];
    embed_X(:,:,n) = embed_Xn;
    
end

%%%Step 3:Construct local constant coefficient SDE simulator at net point
%%%		dXt = b_n dXt + sigma_n dW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T,c,b,sigma,mu] = construct_local_SDE(net_info,embed_L,embed_X,p,t0,m,d);

%%%Save new simulator parameters for ATLAS
new_S.T = T;
new_S.b = b;
new_S.c = c;
new_S.sigma = sigma;
new_S.mu = mu;
new_S.Phi = Phi;

save('current_atlas.mat','new_S','L','embed_L');

if testing
    %%%Save data collected on simulation distances and new t0 values
    atlas_diagnostics.mean_net_distances = mean_net_distances;
    save('current_atlas.mat','atlas_diagnostics','-append');
    
    if vary_t0
        save('current_atlas.mat','t0_list','-append');
    end
    
end

if timestamp_save
    filename =[datestr(now, 'dd_mmm_yyyy_HH_MM'),'_','d',num2str(d),'_','construction'];
    save(filename);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,c,b,sigma,mu] = construct_local_SDE(net_info,embed_L,embed_X,p,t0,m,d)
%TODO: change comments, currently for an older function
%inputs:
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
%outputs:
%		T - d X d X N X N array containing transition mappings between charts,
%		    T(:,:,i,j) is a d x d matrix which shifts chart coords from chart i to j
%		c - d X N X N array containing local coordinates of neighboring charts,
%			 c(:,i,j) is a d x 1 vector, isj-th neighbor of net point n,
%		     expressed in chart n coords(LMDS wrt n's landmarks)
%		b - d X 1 vector, drift coordinate for n's local SDE
%		sigma - d x d x N array, sigma(:,:,n) is d x d array, diffusion coords
%		for n's local SDE
%		mu - d X N X N array containing mean landmarks of local neighboring charts,
%		    mu(:,i,j) is a d x 1 vector, the average of chart j's landmarks
%		    expressed in chart n coords

%Note: following notation in ATLAS paper, 'y_n' is the embedded net point n

%%%Read in net_info and chart t embeddings
%Note: embed_Ln = [y_n |y_n landmarks |nbr1 |nbr1 landmarks |nbr2|nbr2 landmarks...];

neighbors = net_info.neighbors;
deg = net_info.deg;
max_deg = net_info.max_deg;
N = size(net_info.net,2);

T = zeros(d,d,N,N); 		%T is the collection of transition maps from chart to chart
c = zeros(d,N,N); 			%c is the collection of chart centers for each chart
b = zeros(d,N); 			% collection of drift coeff. for local SDEs
sigma = zeros(d,d,N); 		% collection of diffusion matrices for local SDEs;
mu = zeros(d,N,N); 			%mu is the collection of average landmark for each chart

parfor n = 1:N
    fprintf("Computing local SDE for point %d of %d\n",n,N);
    num_nbr_self = deg(n) + 1;
    net_nbr = [n neighbors(n,1:deg(n))];
        
    embed_Ln = embed_L(:,1:(m+1)*(num_nbr_self),n);
    embed_Xn = embed_X(:,:,n);
    
    %%%Compute chart c_n for y_n and center around embedded y_n
    %NOTE: indexing to find nbr is based on that
    %net_nbr = [n nbr1 nbr2 ]
    %embed_Ln = [y_n landmarks| y_nbr1 landmarks |y_nbr2 ...]
    centers_row = zeros(d,N);
    centers_row(:,net_nbr) = embed_Ln(:,1:m+1:(num_nbr_self -1)*(m+1) + 1);
    c(:,n,:) = centers_row;
       
    %%%Compute diffusion, drift coefficients around y_n
    b(:,n) = (1/ (p*t0) )*sum(embed_Xn,2);      %drift for y_n
    sigma(:,:,n) = sqrt(1/t0)*sqrtm(cov((embed_Xn).')); %diffusion for y_n
    
    %%%Compute switching maps
    %NOTE: indexing to find nbr is based on that
    %net_nbr = [n nbr1 nbr2 ]
    %embed_Ln = [y_n landmarks| y_nbr1 landmarks |y_nbr2 ...]
    
    %For parallelizing, storing all updates for the nbr loop in arrays
    mu_nj_row = zeros(d,N);
    T_nj_row = zeros(d,d,N);
    
    for i = 1 : num_nbr_self
        j = net_nbr(i);	%global index of neighbor of n
        
        %%%Find embeddings of chart overlaps between chart for net point n and nbr net point j
        %L_nj is Phi_n(A_n U A_j), the embedding of landmarks for net points n,j in chart n's coordinates
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
                
        % Collect landmark mean, transition map of neighbor for storage
        mu_nj_row(:,j) = mu_nj;
        T_nj_row(:,:,j) = T_nj;
        
    end
    % end transition map computing for neighbor
    
    
    % Collect landmark mean, transition map of all neighbors for storage
    mu(:,n,:) = mu_nj_row;
    T(:,:,n,:) = T_nj_row;
    
end
% end transition map computing for neighbor

end
%% Testing and Diagnostic Functions

function test_landmark_index(n,m,N,neighbors,deg,A,Ln)
%Checking that landmark data struct Ln is properly organized
%inputs:
%		n - current net point, function constructs local_SDE for chart n
%       m - number of landmarks per delta-net point
%       N - number of points in delta net
%		neighbors - N x max_deg array, row j is neighbor indices of net	point j
%		deg - N x 1 vector, entry j is number of neighbors(degree) of net point j
%		embed_L - LMDS applied to L_n, d x p array
%		embed_X - LMDS applied to X_n, d x num_nbr array
%       Ln - D X num_nbr*m array, columns are the landmarks associated with delta net point net(:,n)
%           as well as the neighbors of net(:,n)

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

function [avg_dist_to_net] = diagn_mean_net_dist(n,D,net,delta,t0,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diagnostic function to check average distances of sample paths from an
%   invidiual net point, and compare with delta value
%Ideally, E[X - net(:,n)] ~ delta
%If this isn't the case, consider changing the sim time t0
%
%inputs:
%		n - current net point, function constructs local_SDE for chart n
%		net - D X N array, columns are data points of delta net
%		t0 - desired time of simulation for each call to S
%       X - D x p array, sample paths from net point(:,n)
%               column i is endpoint of a simulation of time t0
%               starting at delta-net point net(:,n)
%outputs:
%		avg_dist_to_net: approx of E[X - net(:,n)]

%TODO: replace euclidean norm with more general norm

if D == 1
    avg_dist_to_net =  mean(abs(X(:,:,n) - net(:,n)),2);
else
    avg_dist_to_net =  mean(vecnorm(X(:,:,n) - net(:,n)),2);
end
fprintf("avg distance from X to net point n:	%f \n",avg_dist_to_net);
fprintf("	delta: 	%f\n",delta);
fprintf("	t_0: 	%f\n",t0);

end
