pot = @(x)CMpot2D(x);
fun = @(x)(-CMgrad2D(x)); % dx = -grad(U(x))dt + dW
delta = 0.2;

MakeDeltaNet(delta,pot,fun);


