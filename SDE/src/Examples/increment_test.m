function [times] = increment_test(data_size,dt_sim)

increment = floor(data_size/1e6)
times = dt_sim*[1e6*ones(1,increment),mod(data_size,1e6)]
for k = 1:length(times)
	times(k)
end

end
