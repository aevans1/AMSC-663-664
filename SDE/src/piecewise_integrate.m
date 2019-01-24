function piecewise_integrate(dx,init_x,init_y,d)
%Integrates a piecewise-linear function from list of differences and plots
%inputs: dx - stepsize in x direction
%        init_x - first x-coord corresponding to first differences val
%        init_y - same as above, for y
%        d - array [d_1 d_2 ... ] of difference values
%outputs: plot of piece-wise linear function constructed with inputs

%%%Initialize variables
n = size(d,2);
x = zeros(1,n);
y = zeros(1,n);

%Start in the middle of x-domain
mid = round(n/2);

x(mid) = init_x;
y(mid) = init_y;

i = 0;
while mid - i > 1 && mid + i < n
%%%%Continue integrating left and right from middle until edge reached

	%step to right...
	x(mid+ i + 1) = dx + x(mid + i);
	y(mid+ i + 1) = y(mid + i) + d(mid + i)*dx;

	%step to the left..
	x(mid - i - 1) = x(mid - i) - dx;
	y(mid - i - 1) = y(mid - i) - d(mid - i)*dx;
    
    i = i + 1;
end

plot(x,y);

end


%Old version of the plotter, starts at left-hand of boundary, moves to the
%right

% function simple_plot(dx,init_x,init_y,d)
% x = [];
% y = [];
% x(1) = init_x;
% y(1) = init_y;
% for i = 1:size(d,2)
% 	x(i+1) = dx + x(i);
% 	y(i+1) = y(i) + d(i)*dx;
% end
% plot(x,y);
% 
% end

