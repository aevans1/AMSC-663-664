n = 500;
for i = 1:n
	ans = norm(rand()*(new_Lk-mu(1,9).mu) - (new_Lj-mu(9,1).mu))
	if (ans < norm(T(1,9).T*(new_Lk-mu(1,9).mu) - (new_Lj-mu(9,1).mu)))
		fprintf("lower!");
		break
	end
end

