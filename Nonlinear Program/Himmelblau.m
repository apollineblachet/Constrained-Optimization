clear all 
close all
clc

%% Plot Himmelblau problem

x = -5:0.005:5;
y = -5:0.005:5;
[X,Y] = meshgrid(x,y);
F = (X.^2+Y-11).^2 + (X + Y.^2 - 7).^2;
v = [0:2:10 10:10:100 100:20:200];
[c,h]=contour(X,Y,F,v,'linewidth',2);
colorbar
yc1 = (x+2).^2;
yc2 = (4*x)/10;
hold on
    fill(x,yc1,[0.7 0.7 0.7],'facealpha',0.2)
    fill([x x(end) x(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.2)
hold off
xlim([-5,5])
ylim([-5,5])