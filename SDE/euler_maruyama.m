%Euler-Maruyama method for linear SDE
%SDE is of form dx = f(x) dt + mu*X dW, X(0) = Xzero
%Currently using 1-D two-well potential U, f(x) = -grad(U)
%Parameters:
	%mu: controls random portion of sde above
	%Xzero: initial condition for SDE, in [0,1], e.g. 0.3
	%T: right endpoint of time domain [0,T], e.g. 100
	%N: number of timesteps, e.g: 2^15
	%R: positive integer, ratio of X timesteps to W timesteps, e.g. 4
	%det: boolean, 1 means compute deterministic DE, 0 means don't

det = 1;

%%%Compute Brownian path
rng(sum(clock));
mu = 1; Xzero = 0.3;
T = 100; N = 2^15; dt = T/N;
dW = sqrt(dt)*randn(1,N);
W = cumsum(dW);

%%%Euler-Maryama parameters for SDE
R = 1; Dt = R*dt; L = N/R;
Xem = zeros(1,L);
Xtemp = Xzero;

%%%Euler parameters for  DE dx/t = f(x)
if det
	X = Xzero;
	Xdet = zeros(1,L);
end

%%%Compute SDE(optional: compute DE alongside)
for j = 1:L
	Winc = sum(dW(R*(j-1)+1:R*j));	%increment brownian motion
	f = -32*Xtemp*(Xtemp-1)*(2*Xtemp-1);
	Xtemp = Xtemp + Dt*f + mu*Winc;
	Xem(j) = Xtemp;

	%%%Compute deterministic DE dx/dt = f(x) 
	if det
		f_det = -32*X*(X-1)*(2*X-1);
		X = X + Dt*f_det;
		Xdet(j) = X;
	end

end


%%%Plotting
plot([0:Dt:T],[Xzero,Xem],'b'), hold on
xlabel('t');
ylabel('X');

if det
	plot([0:Dt:T],[Xzero,Xdet],'r');
end

%emerr = abs(Xem(end) - Xtrue(end))

