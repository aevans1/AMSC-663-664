%%%Parameters and Equations take from Lv et al 2014 Paper
%%%Simulation in 3D of perturbed ODE for the mean field ODE satisfied by mRNA,
%proteins and dimers

dt = 0.0001;
T = 10; 
N = floor(T/dt); %num timesteps
D = 3; %dimension of problem
epsilon = 1; %noise parameter


%xi - asympt. stable equilibria representing inactive switch
%xa - asympt. stable equilibria representing active switch
%xs - morse index one saddle

xi = [0.0402067142317042; 1.60826856926817; 0.000258652779089588],
xa = [29.3768600805981; 1175.07440322392; 138.079985311206],
xs = [10.5829; 423.3173; 17.9198]


X = zeros(D,N+1); 

Xzero = xi;
%Xzero = [0;0;0];
%Xzero = rand(3,1);

X(:,1) = Xzero;

%%%Compute SDE until time dt*N
Xtemp = Xzero;
for j = 1:N
	dW = sqrt(dt)*randn(D,1); %Brownian random increment in R^D
	X(:,j+1) = X(:,j) + switchgrad(X(:,j))*dt + sqrt(epsilon)*dW;
end

%Xtrunc = X(:,1:1000:end);

plot3(X(1,:),X(2,:),X(3,:));
%plot3(Xtrunc(1,:),Xtrunc(2,:),Xtrunc(3,:));
hold on;




hold on;
plot3(xi(1),xi(2),xi(3),'.','Color','r','MarkerSize',10);
%plot3(xa(1),xa(2),xa(3),'.','Color','g','MarkerSize',10);
%plot3(xs(1),xs(2),xs(3),'.','Color','y','MarkerSize',10);





function dX = switchgrad(x)
%%%Parameters
k0 = 1; %DNA activation
k1 = 0.0002; %rate of dimerization
gamma_1 = 2;  %rate of de-dimizeration
gamma_0 = 50; %DNA inactivation
gamma_m = 10; %mRNA decay
gamma_n = 1; %protein decay
a = 400;	%transcription
a0 = 0.4;
b = 40;		%translation


%x is a vector in R^3
m = x(1);
n = x(2);
d = x(3);

%
dm = (a0*gamma_0 + a*k0*d)/(gamma_0 + k0*d) - gamma_m*m;
dn = b*m - gamma_n*n - 2*k1*n^2 + 2*gamma_1*d;
dd = k1*n^2 - gamma_1*d;

dX = [dm;dn;dd];

end
