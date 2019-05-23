function makestats_tpaths_gswitch()
%Load up transition data for genetic switch and make bar graphs
%comparing ATLAS transition data with original simulator transition data

%% Read in Data

%%%read in switch times from ATLAS
data = load('learned_tswitch.mat');
learned_switches = data.learned_switches;

%%% read in switch times from original simulator
data = load('original_tswitch.mat');
original_switches = data.original_switches;

num_regions = 2;

%% Plot Bar Graphs

fprintf('ATLAS:\n');
Tbar = MeanSwitchTimes(num_regions,learned_switches);
fprintf('Original\n');
Tbar0 = MeanSwitchTimes(num_regions,original_switches);

figure(2); clf; hold on;
bar([Tbar,Tbar0]);
xticks([1:6]);
xticklabels({'1->2','2->1'});
set(gca,'Fontsize',20);
ylabel('Mean transition time','Fontsize',20);
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
end


