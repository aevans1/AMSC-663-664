for n = 1:size(net,2)
	sim_net(:,n) = S(net(:,n),1,0.1);
end
save('sim_net');
