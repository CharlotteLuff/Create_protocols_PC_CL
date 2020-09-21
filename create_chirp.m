close all
clear all
clc

dt = 0.004;
f0 = 0.5;
t1 = 10;
f1 = 10;
t =  -t1+2:dt:t1;

I1 = chirp(t,f0,t1,f1,'quadratic',[],'concave');
figure, plot(I1)
saveas(gcf,['chirp_I1_2ep' '.fig']);

fx0 = 0.5;
tx1 = 10;
fx1 = 10;
tx =  -tx1+2:dt:tx1;

I2 = chirp(tx,fx0,tx1,fx1,'quadratic',[],'convex');
figure, plot(I2)
saveas(gcf,['chirp_I2_2ep' '.fig']);

figure, plot(I1 + I2)
save('chirp_2ep' ,'I1','I2')
saveas(gcf,['chirp_FM_2ep' '.fig']);

I1 = I1 + I2;
I2 = zeros(1,length(I1));
saveas(gcf,['chirp_I2_2ep' '.fig']);

figure, plot(I1)
save('chirp_1ep' ,'I1','I2')
saveas(gcf,['chirp_FM_1ep' '.fig']);