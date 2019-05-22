function calculate_switches_test()
%Two simple hand-solved examples for the switch_times testing,
%tests for function 'calculate_switches', included in this file and used in
%transition_data.m in 'Examples' folder


%% Example 1
%paths: sequences of regions (1 or 2 in this case) visited over time,
%dt: timestep of underlying dynamics, in this case we just set dt = 1
paths = [1 1 1 0 1 0 0 2 2 2 1 0 2 0 1 1 2]; 
dt = 1;
num_steps = 17;

%Organizing transition data:
%ith row: information on ith transition observed
%1st col: initial region of transition
%2nd col: final region of transiton
%3rd col: how many steps for transition
%For 'paths' above
%     1 2 7 -> transition from region 1 to region 2 in 7 steps
%     2 1 3 -> ....
%     1 2 2 
%     2 1 2
%     1 2 2
true_times = [1 2 7; 2 1 3; 1 2 2; 2 1 2; 1 2 2];
if (calculate_switches(paths,dt,num_steps) == true_times)
	fprintf("test 1 correct! \n");
else
	fprintf("test 1 error! \n");
end


%% Example path 2, 
paths = [4 4 4 0 1 1 2 3 1 4 3 3 3];
num_steps = size(paths,2);
%%
%compute true timings of above by visual inspection
true_times = [4 1 4; 1 2 2; 2 3 1; 3 1 1; 1 4 1; 4 3 1];

if (calculate_switches(paths,dt,num_steps) == true_times)
 
	fprintf("test 2 correct! \n");
else
	fprintf("test 2 error! \n");
end

function switch_times = calculate_switches(path_regions,dt,num_steps)
%Compute switching times given a list of regions switched over timesteps
%with prescribed timestep dt
%inputs:
%       path_regions - 1 x num_steps array, entry i is the index of a
%               region visited at time dt*i
%       dt - timestep
%       num_steps - length of path_regions
%outputs:
%       switch_times - array, each row has values [old_region new_region
%           t] corresponding indices to the regions switched and the time
%               t the switch took
%
%%%Initialize
switch_times=[];
switch_count = 0; %number of switches
index = path_regions(1); 
flag = sign(index); %flag = 0 means a region hasn't been visited yet
t = 0;

for n = 2:num_steps
    new_index = path_regions(n);
    
    t = t+dt;
    
    %if new region is visited, restart the timer
    %   if we've left region 0, count a switch and update switch list
    if new_index > 0              
        if flag == 0
            index = new_index;
            t = 0;
            flag = 1;
        elseif index ~= new_index
            switch_count = switch_count + 1;
            %fprintf("Transition %d \n",switch_count);
            %we've observed a transition!
            %switch_count = switch_count+1
            switch_times = [switch_times;index,new_index,t];
            t = 0;
            index = new_index;
        end
    end
    
end
end

end