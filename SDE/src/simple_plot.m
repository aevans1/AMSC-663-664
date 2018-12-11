function simple_plot(dx,init_x,init_y,d)
x = [];
y = [];
x(1) = init_x;
y(1) = init_y;
for i = 1:size(d,2)
	x(i+1) = dx + x(i);
	y(i+1) = y(i) + d(i)*dx;
end
plot(x,y);

end
%function simple_plot(x,y)
%
%end
