function U = CMpot2D(x)
% 2D potential from Crosskey&Maggioni
p1 = [0;0];
p2 = [1.5;0];
p3 = [0.8;1.05];
p = [p1,p2,p3];
c = [0.2,0.2,1/6];
for i = 1 : 3
   d(:,i) = x - p(:,i);
   r(i) = d(:,i)'*d(:,i)/c(i);
   e(i) = exp(-r(i));
end
% potential
U = -log(sum(e));
end
