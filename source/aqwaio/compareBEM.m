clc; clear all; close all;

c = load('./aqwa.mat');
c = c.hydroData(1); % take first body only
c.fExtRe = c.fExt.Re;
c.fExtMag = c.fExt.Mag;
c.fExtIm = c.fExt.Im;

m = load('./aqwa-data-wecSimHydroData1.mat');

i = 1;
plot(c.period,c.fExtMag(i,:),'o-'); hold on
plot(m.period,m.fExtMag(i,:),'x-');
legend('carlos','mike')

figure;
i = 1;
plot(c.period,c.fExtRe(i,:),'o-'); hold on
plot(m.period,m.fExtRe(i,:),'x-');
legend('carlos','mike')

figure;
i = 1;
plot(c.period,c.fExtIm(i,:),'o-'); hold on
plot(m.period,m.fExtIm(i,:),'x-');
legend('carlos','mike')

figure;
plot(c.period); hold on;
plot(m.period')