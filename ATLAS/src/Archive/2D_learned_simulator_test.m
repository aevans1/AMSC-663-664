%NOTE: this is the main test file for comparing simulators
%Data generated by this file should be used by a number of other test files,
%one for transition paths, one for histograms...etc

setting = 'transition_paths';
%%%TESTING: binning for original simulator
if ~exist('setting','var'), setting = 'sim_compare'; end
%%

setting_str = sprintf("%s",setting);
filename =[datestr(now, 'dd_mmm_yyyy_HH_MM'),'_','d',num2str(d),'_','learned_simulator_test','_',setting];
fileID = fopen([filename,'.txt'],'w');

fprintf("Temporary: setting seed \n");
seed = 1;
fprintf("Seed: %d \n",seed);

rng(seed);
p = 1;
T = 500;

%TESTING
fprintf("Testing transition path tracking \n");
%NOTE: both of these parameters are very problem dependent
%TODO: move these params to another file?
regions = [0 1.5 0.8; 0 0 1.05];
dist = (1/4);

%Two ways: random sample initial locations, or all net points
%1) Random samples
num_locations = 12;

%2) Sample from net
%num_locations = size(net,2);

%initialize sample collection and transition collections
inverse_collection = zeros(d,num_locations*p);
X_collection = zeros(d,num_locations*p);
Y_switch_collection(1,1).list = [];
X_switch_collection(1,1).list = [];
for i = 1:length(regions)
	for j = i + 1: length(regions)
		if j~=i
			%append new i,j switches to existing list	
			Y_switch_collection(i,j).list = [];
			Y_switch_collection(j,i).list = [];
			X_switch_collection(i,j).list = [];
			X_switch_collection(j,i).list = [];
		end
	end
end


fprintf(fileID,'example: %d\n',example);
fprintf(fileID,'\n');
fprintf(fileID,'delta=%f\n',delta);
fprintf(fileID,'m=%d\n',m);
fprintf(fileID,'t_0=%f\n',t_0);
fprintf(fileID,'dt(ATLAS)=%f\n',dt);
fprintf(fileID,'dt(original)=%f\n',dt_original);
fprintf(fileID,'p=%d\n',p);
fprintf(fileID,'T=%f\n',T);
fprintf(fileID,'\n');

fprintf(fileID,"Parameters from simulator testing: \n");
fprintf(fileID,'regions = %d\n',regions);
fprintf(fileID,'dist=%f\n',dist);
fprintf(fileID,'num_locations=%d\n',num_locations);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TESTING: where to set Xzero?
%NOTE: doesn't seem to affect transition results
%Xzero = 1.6*rand(d,1) - 0.3;
for m = 1:num_locations
	
%	1)Uniform sample from domain:
%
%	Example 1 or 2: uniform from [-0.3,1.3]
	%Xzero = 1.6*rand(d,1) - 0.3;

	%Example 3: uniform from [-1,2.5] x [-1, 2]
	Xzero(1,1) = 3.5*rand() - 1;
	Xzero(2,1) = 3*rand() - 1;

% 	2)Use net point
	%Xzero = net(:,m);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Reduced(ATLAS) Simulator
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%Simulate Y paths
%	Find the initial start to start reduced simulator
	for n = 1:size(net,2)
	   	distances(n) = norm(Xzero - net(:,n));
	end
	init_chart = find(distances == min(distances),1);

% 	Embed initial condition for original simulator so it can berun in reduced
% 	simulator
	[local_L,Yzero] = LMDS(L(init_chart).L,Xzero,rho,d);
	
% 	Center projected point around the chart center
	Yzero = Yzero - local_L(:,1);
	
	  [Ypaths,charts,Y_regions_visited] =learned_simulator(Yzero,p,dt,T,new_S,neighbors,d,delta,net,init_chart,regions,dist);

 	%%%%Map reduced sim paths back to delta net in original space
	inverse_Ypaths = [];
	for i = 1:p
		inverse_Ypaths(:,i) = net(:,charts(i));
	end
	inverse_collection(:,1 + (m-1)*p: m*p) = inverse_Ypaths;

	%%%Collect transition data between regions of interest
	[Y_switch_times] = transition_data(regions,Y_regions_visited);
	for i = 1:length(regions)
		for j = i + 1: length(regions)
			if j~=i
				%append new i,j switches to existing list	
				switch_list = Y_switch_collection(i,j).list;
				new_switches = Y_switch_times(i,j).list;
				switch_list = [switch_list new_switches];
				Y_switch_collection(i,j).list = switch_list;

				%append new j,i switches to existing list	
				switch_list = Y_switch_collection(j,i).list;
				new_switches = Y_switch_times(j,i).list;
				switch_list = [switch_list new_switches];
				Y_switch_collection(j,i).list = switch_list;
			end
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Original Simulator
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%Simulate X paths
	init_x = Xzero;
	[X,X_regions_visited] = simulator(Xzero,m,T,f,dt_original,regions,dist);
	
	%%%bin the endpoints into delta net, collect data
	X_delta = [];
	for i = 1:p
		for n = 1:size(net,2)
	    	distances(n) = norm(X(:,i) - net(:,n));
		end
		j =find(distances == min(distances),1);
		X_delta(:,i) = net(:,j);
	end
	X_collection(:,1 + (m-1)*p: m*p) = X_delta;	

	%%%Collect transition data between regions of interest
	[X_switch_times] = transition_data(regions,X_regions_visited);
	for i = 1:length(regions)
		for j = i + 1: length(regions)
			if j~=i
				%append new i,j switches to existing list	
				switch_list = X_switch_collection(i,j).list;
				new_switches = X_switch_times(i,j).list;
				switch_list = [switch_list new_switches];
				X_switch_collection(i,j).list = switch_list;

				%append new j,i switches to existing list	
				switch_list = X_switch_collection(j,i).list;
				new_switches = X_switch_times(j,i).list;
				switch_list = [switch_list new_switches];
				X_switch_collection(j,i).list = switch_list;
			end
		end
	end
	fprintf("Finished samples for point %d of %d \n",m,num_locations);
end

%%%Collect all average transitions for simulators
for i = 1:length(regions)
	for j = i + 1: length(regions)
		if j~=i
			X_switch_averages(i,j) = mean(X_switch_collection(i,j).list);
			X_switch_averages(j,i) = mean(X_switch_collection(j,i).list);

			Y_switch_averages(i,j) = mean(Y_switch_collection(i,j).list);
			Y_switch_averages(j,i) = mean(Y_switch_collection(j,i).list);
		end
	end
end


switch setting
	case 'sim_compare'
		edges = [-0.5:delta:1.5];
		
		figure;
		mu_hat = histogram(inverse_collection,edges);
		hold on;
		
		mu = histogram(X_collection,edges);
		legend('reduced_sim','original_sim');
	
	case 'transition_paths'
	%applies for either of example 1 or 2
	%TODO: print statement with example number, relevant parameters, etc.
	
	%TODO: more robust outputs, bar graphs, etc.	
	%fprintf("Original Sim, avg transition 1-->2: %f\n",X_switch_averages(1,2));
	%fprintf("Original Sim, avg transition 2-->1: %f\n",X_switch_averages(2,1));
	%fprintf("ATLAS Sim, avg transition 1-->2: %f\n",Y_switch_averages(1,2));
	%fprintf("ATLAS Sim, avg transition 2-->1: %f\n",Y_switch_averages(2,1));

	X_switch_bar = X_switch_averages.';
	X_switch_bar = X_switch_bar(:);
	Y_switch_bar = Y_switch_averages.';
	Y_switch_bar = Y_switch_bar(:);
	figure;
	hold on;
	bar([X_switch_bar,Y_switch_bar]);
	xticks([1:6]);
	xticklabels({'1->2','1->3','2->1','2->3','3->1','3->2'});
	set(gca,'Fontsize',20);
	ylabel('Mean transition time','Fontsize',20);

end
save(filename);
fclose(fileID);

