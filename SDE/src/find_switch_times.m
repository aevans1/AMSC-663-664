function find_switch_times()
% Runs the original simulator and finds mean switch times between wells
% in Crosskey-Maggioni 3-well potential
global p1 p2 p3 ppp wrad
p1 = [0;0];
p2 = [1.5;0];
p3 = [0.8;1.05];
ppp = [p1,p2,p3];
%%
pot = @(x)CMpot2D(x);
fun = @(x)(-CMgrad2D(x)); % dx = -grad(U(x))dt + dW
data = load('delta_net_CM2D.mat');
net = data.net; % coordinates of vertices of the delta-net
E = data.E; % edges of the delta-net
xgrid = data.xg;
ygrid = data.yg;
Ugrid = data.U;
deg = data.deg;
maxdeg = data.maxdeg;
nei = data.nei;
% %% graphics
% figure(3); clf; hold on;
% contour(xgrid,ygrid,Ugrid,linspace(0,3.7,20),'Linewidth',1);
% xlabel('x','Fontsize',20);
% ylabel('y','Fontsize',20);
% set(gca,'Fontsize',20);
% drawnow;
%
data = load('TswitchCM2D.mat');
tmax = data.Nsteps*data.dt;
fprintf('The total simulation time is %d\n',tmax);
wrad = data.wrad;
%%
dt = 0.005; % timestep for Euler-Maruyama
%% dx = fun(x)*dt + dw
nsteps = round(tmax/dt);
x = net(:,1);
iw  = well(x);
flag = sign(iw);
dw = sqrt(dt)*randn(2,nsteps);
Tswitch = [];
t = 0;
for i = 1 : nsteps
    x = x + fun(x)*dt + dw(:,i);
    t = t + dt;
    iwnew = well(x);
%     fprintf('well = %d\n',iwnew);
%     plot(x(1),x(2),'.','Markersize',10); drawnow;
    if iwnew > 0 % got into some well
        if flag == 0 % this is the first well
            iw = iwnew;
            flag = 1;
            t = 0;
        else
            if iw ~= iwnew
                Tswitch = [Tswitch; [iw,iwnew,t]];
                iw = iwnew;
                t = 0;
            end
        end
    end
end
save('OriginalTswitch.mat','Tswitch');
end

%%
function iw = well(x)
global p1 p2 p3 ppp wrad
[dmin,imin] = min([norm(x - p1),norm(x - p2),norm(x - p3)]);
if norm(x - ppp(:,imin)) < wrad
    iw = imin;
else
    iw = 0;
end
end




