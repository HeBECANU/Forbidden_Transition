% data = [20.5,156
% 31.9,375
% 11.42,277];
% data = [35.1,870
% 20.8,322
% 0.4,110];
data=[36.4,1100
23.6,740
35.9,1080
7.47,220
0.480,103
22.3,886];
% data=[21.7,705
% 41,1170];

y = data(:,1);
x = data(:,2);

X = [ones(length(x),1) x];
b = X\y;
figure(111);
%clf
scatter(x,y,'kx')
yCalc2 = X*b;
hold on
plot(x,yCalc2,'-')

