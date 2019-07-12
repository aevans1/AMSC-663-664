function [times,averages] = transition_data(regions,paths,dt)
%Given a list of regions visited (paths), e.g. [1 2 0 1 1 2],
%and specified regions of interest,e.g. [1 2], returns transition times between
%regions
%inputs:  paths - row vector of integers
%		  regions - subset of possible entries in paths
%outputs: averages - array, entry (i,j) is mean transition time from region
%		  			 index i to region index j
%		  times - struct array, struct times(i,j).list gives list of all
%		  		  transition times between region index i and region index j in paths         

%NOTE: regions length should be bigger than 1, else nothing gets returned
for i = 1:length(regions)
	for j = i+1:length(regions)
		i_to_j =switch_times(i,j,paths,dt);
		times(i,j).list = i_to_j;
		averages(i,j) = mean(i_to_j);
		%fprintf("transition times for %d to %d:\n",regions(i),regions(j));
		%fprintf("%d ",i_to_j);
		%fprintf("\n");
		%fprintf("average time from %d to %d: %f \n",regions(i),regions(j), mean(i_to_j));
		%fprintf("\n");
		
		j_to_i = switch_times(j,i,paths,dt);
		times(j,i).list = j_to_i;
		averages(j,i) = mean(j_to_i);
		%fprintf("transition times for %d to %d \n",regions(j),regions(i));
		%fprintf("%d ",j_to_i);
		%fprintf("\n");
		%fprintf("average time from %d to %d: %f \n",regions(j),regions(i), mean(j_to_i));
		%fprintf("\n");
	end
end
end
function times = switch_times(a,b,paths,dt)
%Given a list of regions visited (paths), e.g. [1 2 0 1 1 2],
% and regions a, b e.g. a = 1, b = 2, returns list of transition times from a
% to b in paths
%inputs:  paths - row vector of integers
%		  a - starting region for transition, should correspond to an entry in
%		  	  paths
%		  b - ending region for transition, should correspond to an entry in
%		  	  paths
%		  times - row vector of all transition times between region index i
%		  		  and region index j in paths         
	i = 1;
	j = 1;
	times = [];

	%%%loop through all possible left endpoints a, check for right endpoints b
	while i < length(paths)

		%iterate until left endpoint is a or out of bounds
		while paths(i)~= a && i < length(paths)
			i = i + 1;
		end
	
		%choose right endpoint for transition
		j = i + 1;
		
		if j <=length(paths)
			%iterate until right endpoint is b or out of bounds
			while paths(j) ~= b && j < length(paths)
				j = j + 1;
			end
			
			%if no out of bounds from earlier, add time difference to list
			if (paths(i) == a) && (paths(j) == b)
				times(end + 1) = j - i;
			end
			i = j + 1;	
		else
			i = length(paths); %break
		end
    end
    times = dt*times; %normalize to get actual times
end
