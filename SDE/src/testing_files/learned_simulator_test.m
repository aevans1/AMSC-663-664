%%%%TESTING: binning for original simulator
%fprintf("temporarily reducing number of calls to learned simulator \n");

%TEMPORARY
fprintf("Temporary: setting seed");
rng(1);

p = 1000;

Xzero = rand(d,1);
for n = 1:size(net,2)
   	distances(n) = norm(Xzero - net(:,n));
end
init_chart = find(distances == min(distances),1);

%Embed initial condition for original simulator so it can berun in reduced
%simulator
[~,Yzero] = LMDS(L(init_chart).L,Xzero,rho,d);

%Center projected point around the chart center
Linit = L(init_chart).L;
Yzero = Yzero - Linit(:,1);

T = 0.2;

[Ypaths,charts] =learned_simulator(Yzero,p,dt,T,new_S,neighbors,d,delta,net,init_chart);

inverse_Ypaths = [];
for i = 1:p
	inverse_Ypaths(:,i) = net(:,charts(i));
end

figure;
histogram(inverse_Ypaths,10);
hold on;

init_x = Yzero;
%simulate p paths around net point y_n
X = S(init_x,p,T);

%%bin the endpoints into delta net
X_delta = [];
for i = 1:p
	
	for n = 1:size(net,2)
    	distances(n) = norm(X(:,i) - net(:,n));
	end
	j =find(distances == min(distances),1);
	X_delta(:,i) = net(:,j);
end

histogram(X_delta,10);



