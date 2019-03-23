function a = simulator(y,m,tmax,fun)
dt = 0.005; % timestep for Euler-Maruyama
%% dx = fun(x)*dt + dw
nsteps = round(tmax/dt);
for j = 1 : m
    dw = sqrt(dt)*randn(2,nsteps);
    x = y;
    for i = 1 : nsteps
        x = x + fun(x)*dt + dw(:,i);
    end
    a(:,j) = x;
%     if mod(j,100) == 0
%         fprintf('simulator: j = %d\n',j);
%     end
end
end