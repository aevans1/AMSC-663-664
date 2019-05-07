%Three simple hand-solved examples for the switch_times testing
%true_times is the ground truth for how long each switch time is between
%regions

paths = [1 1 1 0 1 0 0 2 2 2 1 0 2 0 1 1 2];
regions = [1 2];
true_times = [];

true_times(1,2).list = [7 2 2];
true_times(2,1).list = [3 2];
[times,averages] = transition_data(regions,paths)
tf = isequaln(times,true_times)
if tf
	fprintf("Correct! \n");
else
	fprintf("Error! \n");
end

paths = [4 4 4 0 1 1 2 3 1 4 3 3 3];
regions = [1 2 3 4];
true_times = [];
true_times(1,2).list = [2];
true_times(2,1).list = [2];
true_times(1,3).list = [3 2];
true_times(3,1).list = [1];
true_times(1,4).list = [5];
true_times(4,1).list = [4];
true_times(2,3).list = [1];
true_times(3,2).list = [];
true_times(2,4).list = [3];
true_times(4,2).list = [6];
true_times(3,4).list = [2];
true_times(4,3).list = [7 1];
[times,averages] = transition_data(regions,paths)
tf = isequaln(times,true_times);
if tf
	fprintf("Correct! \n");
else
	fprintf("Error! \n");
end
paths = repmat([1 2 3],1,100);
regions = [1 2 3];
true_times = [];
true_times(1,2).list = repmat([1],1,100)
true_times(2,1).list = repmat([2],1,99)
true_times(1,3).list = repmat([2],1,100);
true_times(3,1).list = repmat([1],1,99);
true_times(2,3).list = repmat([1],1,100);
true_times(3,2).list = repmat([2],1,99);
[times,averages] = transition_data(regions,paths)
tf = isequaln(times,true_times);
if tf
	fprintf("Correct! \n");
else
	fprintf("Error! \n");
end
