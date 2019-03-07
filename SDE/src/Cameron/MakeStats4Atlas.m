function MakeStats4Atlas()
dat = load('TswitchCM2D.mat');
Tswitch = dat.Tswitch;
Tswitch(1,:) = [];
ch = dat.ch;
% Tswitch is the list of switch times between wells
% ch is the sequence of visited charts
%% compare with the original simulator
dat = load('OriginalTswitch.mat');
Tsw0 = dat.Tswitch;
%%
Nwell = 3;
fprintf('ATLAS:\n');
Tbar = MeanSwitchTimes(Nwell,Tswitch);
fprintf('Original\n');
Tbar0 = MeanSwitchTimes(Nwell,Tsw0);


figure(2); clf; hold on;
bar([Tbar,Tbar0]);
xticks([1:6]);
xticklabels({'1->2','1->3','2->1','2->3','3->1','3->2'});
set(gca,'Fontsize',20);
ylabel('Mean transition time','Fontsize',20);
end

%% find mean switch times
function Tbar = MeanSwitchTimes(Nwell,Tswitch)
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


