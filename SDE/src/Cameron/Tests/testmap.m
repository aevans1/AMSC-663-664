function testmap()
close all
data = load('LearnedSimulator_CM2D.mat');
Lnew = data.Lnew;
T = data.T;
mu = data.mu;
 
L1new = squeeze(Lnew(:,:,1))
L7new = squeeze(Lnew(:,:,7))

figure(4); hold on; grid;
col1 = 'b';
col2 = 'r';
plot(L1new(1,1:11),L1new(2,1:11),'.','color',col1,'Markersize',20);
plot(L1new(1,1),L1new(2,1),'.','color',col1,'Markersize',40);
plot(L1new(1,12:22),L1new(2,12:22),'.','color',col2,'Markersize',20);
plot(L1new(1,12),L1new(2,12),'.','color',col2,'Markersize',40);
daspect([1,1,1]);

figure(5); hold on; grid;
col1 = 'r';
col2 = 'b';
plot(L7new(1,1:11),L7new(2,1:11),'.','color',col1,'Markersize',20);
plot(L7new(1,1),L7new(2,1),'.','color',col1,'Markersize',40);
plot(L7new(1,78:88),L7new(2,78:88),'.','color',col2,'Markersize',20);
plot(L7new(1,78),L7new(2,78),'.','color',col2,'Markersize',40);
daspect([1,1,1]);

ch1pt1 = L1new(:,1:11)
ch1pt2 = L1new(:,12:22)
ch2pt1 = L7new(:,78:88)
ch2pt2 = L7new(:,1:11)

T21 = squeeze(T(:,:,7,1))
T12 = squeeze(T(:,:,1,7))
mu21 = squeeze(mu(:,7,1));
mu12 = squeeze(mu(:,1,7));

map12pt1 = ((ch1pt1 - mu12)'*T12)' + mu21;
figure(5);
plot(map12pt1(1,:),map12pt1(2,:),'.','color','k','Markersize',30);

end