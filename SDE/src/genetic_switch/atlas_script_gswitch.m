data = load('reduced_QPotLevelSets.mat');
reduced_patches = data.reduced_patches;

%reduced_patches:
%f/v 1 - 4: quasipotential level sets for xi(lower) equilibria
%f/v 5 - 8: quasipotential level sets for xa(upper) equilibria
%f/v 9 - 10: quasipotential level sets with respect to upper equilibria, but
%            covers both equilibria and saddle
%f/v 11: same as 9 - 10, but sparse set and unusable
%fields f0,f1,f2,...,f10: faces
%fields v0,v1,v2,...,v10: vertices (points to be used for the delta net)

figure;
hold on; grid;
view(3);
col = [winter(4);spring(4);[1 0 0];[0 1 0];[0 0 1]];

%initialize point set for delta-net
%columns are vectors in R^3
init = [];

for k = 1:10
    
    %read in path, plot
    my_field = strcat('v0_',num2str(k));
    v = reduced_patches.(my_field);
    my_field = strcat('f0_',num2str(k));
    f = reduced_patches.(my_field);
    %p = patch('Vertices',v,'Faces',f,'Facecolor',col(k,:),'Edgecolor',col(k,:));
    %alpha(0.4);
    
    %
    %     %keep only the saddle portion of the 10th patch
    %     if ismember(k,[10])
    %         col = v(:,2);
    %         ind = find(col > 300 & col < 500);
    %         v = v(ind,:);
    %     end
    init = [init,v.'];
end

size(init)
%%
%
% k = 10;
% eps1 = 5.0;
% new_v = [];
% my_field = strcat('v0_',num2str(k));
% v = reduced_patches.(my_field);
% v = v.';
% for i = 1:size(v,2)
%
%     x = v(:,i);
%
%     for j = 1:10
%         %Eulear-Maruyama step
%         dW = sqrt(dt)*randn(3,1); %Brownian random increment in R^D
%         x = x + switch_drift(x)*dt + sqrt(eps1)*dW;
%     end
%     new_v = [new_v x];
%
%
%
% end

%init = [init,new_v];



num_points = size(init,2);

%TESTING
% fprintf("TESTING: making a much smaller initial point set\n");
% size(init)
% init = init(:,1:6:num_points);
% num_points = size(init,2);

%%%plot before running simulation
plot3(init(1,:),init(2,:),init(3,:),'.');
title('before');
view(3);

%%

%run small simulation for each net point to 'relax' set of init points
dt = 0.001;
eps1 = 5.0;
eps2 = 20.0;
dist = 200;
rho = @(x,y) norm(x - y);
new_init = [];


%isolate set of points that are under-sampled (near the saddle)

for i = 1:num_points
    x = init(:,i);
    
    for k = 1:50
        %Eulear-Maruyama step
        dW = sqrt(dt)*randn(3,1); %Brownian random increment in R^D
        x = x + switch_drift(x)*dt + sqrt(eps1)*dW;
    end
    new_init = [new_init x];

    
    
    if x(3) < 50
        num_simulated = 200;
        eps = eps2;
    else
        num_simulated = 50;
        eps = eps1;
    end
    
    for k = 1:num_simulated
        dW = sqrt(dt)*randn(3,1); %Brownian random increment in R^D
        x = x + switch_drift(x)*dt + sqrt(eps)*dW;
        new_init = [new_init x];
    end
    
    %
    %     %simulate more and noisier for points near and below the saddle
%     for m = 1:20
        %         dW = sqrt(dt)*randn(3,1); %Brownian random increment in R^D
        %         x = x + switch_drift(x)*dt + sqrt(eps)*dW;
        %         for k = 1:num_simulated
        %             %Eulear-Maruyama step
        %             dW = sqrt(dt)*randn(3,1); %Brownian random increment in R^D
        %             y = x + switch_drift(x)*dt + sqrt(eps)*dW;
        %             new_init = [new_init y];
        %         end
        %     end
        %
        %
        
        fprintf("%d of %d\n",i,num_points);
    end
    
    figure; grid; hold on;
    plot3(new_init(1,:),new_init(2,:),new_init(3,:),'.');
    title('after');
    view(3);
    
    %%%Delta net trajectories option
    %Euler-Maruyama step
    %         dW = sqrt(dt)*randn(3,1); %Brownian random increment in R^D
    %         x = x + switch_drift(x)*dt + sqrt(epsilon)*dW;
    
    %         %Add point to net if delta away from net
    %         if far(x,new_init,delta,rho)
    %             N = N + 1;
    %             new_init = [new_init x];
    %         end
    %end
    
    
    %%%Use just new points, or new points with initial points
    %point_set = [init new_init];
    point_set = new_init;
    
    ind = randperm(size(point_set,2));
    point_set = point_set(:,ind);
    
    save('point_set_new.mat','point_set');
    %%
    %
    % data = load('point_set_new.mat');
    % point_set = data.point_set;
    
    %%plot init set after running simulation
    % figure;grid;
    % hold on;
    % plot3(point_set(1,:),point_set(2,:),point_set(3,:),'.');
    % title('after');
    % view(3);
    % %
    %%
    data = load('point_set_new.mat');
    point_set = data.point_set;
    % %%Create delta net
    delta = 5.0;
    %delta = 7.5;
    rho = @(x,y) norm(x - y);
    delta_net(point_set,delta,rho);
    %%
    %output of delta_net ^^ goes to curren't delta_net, make a copy in
    %gswitch_delta_net for further use with gswitch stuff
    
    copyfile('current_delta_net.mat','gswitch_delta_net.mat');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function is_far = far(x,net,delta,rho)
    %%%NOTE: if input x is part of net, this will return false
    %determines if point 'x' is at least delta-far from set of points 'net'
    %inputs: x is vector in R^D, net is D X N matrix, N vectors in net
    %output: boolean, true if point is far from set
    
    is_far = true;
    N = size(net,2);
    k = 1;
    
    %Check is x is far away from net
    while (is_far && k <= N)
        if rho(x,net(:,k)) < delta
            is_far = false; %break!
        else
            k = k+1;
        end
    end
    
    end
    
    function dX = switch_drift(x)
    %%%Parameters
    k0 = 1; %DNA activation
    k1 = 0.0002; %rate of dimerization
    gamma_1 = 2;  %rate of de-dimizeration
    gamma_0 = 50; %DNA inactivation
    gamma_m = 10; %mRNA decay
    gamma_n = 1; %protein decay
    a = 400;	%transcription
    a0 = 0.4;
    b = 40;		%translation
    
    
    %x is a vector in R^3
    m = x(1);
    n = x(2);
    d = x(3);
    
    %
    dm = (a0*gamma_0 + a*k0*d)/(gamma_0 + k0*d) - gamma_m*m;
    dn = b*m - gamma_n*n - 2*k1*n^2 + 2*gamma_1*d;
    dd = k1*n^2 - gamma_1*d;
    
    dX = [dm;dn;dd];
    
    end