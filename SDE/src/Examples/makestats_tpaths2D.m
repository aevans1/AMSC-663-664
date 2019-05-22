function makestats_tpaths2D()
%Load up transition data from 'transition_paths.m' and makes bar graphs
%comparing ATLAS transition data with original simulator transition data

%NOTE: this file is written in the
%examples 3 for atlas_driver.m file which has 3 regions

%% Read in Data

%%%read in switch times from ATLAS
data = load('learned_tswitch.mat');
learned_switches = data.learned_switches;

%%% read in switch times from original simulator
data = load('original_tswitch.mat');
original_switches = data.original_switches;

num_regions = 3;
%% Plot Bar Graphs
fprintf('Original\n');
Tbar0 = MeanSwitchTimes(num_regions,original_switches);

fprintf('ATLAS:\n');
Tbar = MeanSwitchTimes(num_regions,learned_switches);

figure(2); clf; hold on;
bar([Tbar,Tbar0]);
xticks([1:6]);
xticklabels({'1->2','1->3','2->1','2->3','3->1','3->2'});
set(gca,'Fontsize',20);
ylabel('Mean transition time','Fontsize',20);
legend('Original','ATLAS');

end

function Tbar = MeanSwitchTimes(Nwell,Tswitch)
%TODO: comments

T = zeros(Nwell); % mean switch times between wells
for i = 1 : Nwell
    for j = 1 : Nwell
        if j ~= i
            ind = find(Tswitch(:,1) == i & Tswitch(:,2) == j);
            if ~isempty(ind)
                T(i,j) = mean(Tswitch(ind,3));
                fprintf('Mean switch time from well %d to well %d is %d, #switch = %d\n',i,j,T(i,j),length(ind));
            end
        end
    end
end
T = T';
Tbar = T(:);
Tbar([1 : Nwell + 1 : end]) = []; %delete diagonal entries
% Now Tbar must be
% 1-->2; 1-->3; 2-->1; 2-->3; 3-->1; 3-->2
end


