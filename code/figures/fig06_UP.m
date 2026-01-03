% fig06_UP.m
% Reproduces Figure 6 (UP model) from the paper.
% ODE simulation of NLUC/FLUC output for the UP circuit, overlapped with data.

clear; close all;

%% Experimental data (NLUC/FLUC)
time=[0 5*60 10*60 15*60 30*60 60*60 120*60 240*60];
mRNA=[0.190216862 0.250747676 0.410231036 0.416895276 0.602128859 0.72489613 0.817693199 0.922975691];
sd2=[0.03891389	0.04077406 0.088993685 0.224306244 0.108054285 0.220286972 0.102963159 0.257725747];

%% Parameters (as used in the paper)
k1 = 0.0301;
k2 = 0.0001; 
k3 = 26.5878; %per second
KM1 = 16.4217; %percent
KM2 = 270.2505; %micro-molar
k4 = 0.0088;
k5 = 0.01;
KM3 = 0.5;
PCO4 = 0.08;
k6 = 0.0001;
k7 = 0.0009;
k8 = 9.6688; % Kcat for PCO1
KM4 = 19.5138; % KM1 o2 for PCO1
KM5 = 286.6942; % KM2 rap for PCO1
k9 = 0.03;
HRPE1 = 0.02;
HRPE2 = 0.005;
HRPE3 = 0.015;

%% Initial conditions under normoxia (O = 21)
mRNA0 = mRNA(1);
O = 21;

RAP0 = k5*mRNA0*KM3/(k4*HRPE1-HRPE1*k5*mRNA0);

f = @(t,a) [-k3*PCO4*O/(KM1+O)*a(1)/(KM2+a(1))+k1-k2*a(1)+k9*a(1)*HRPE3/(KM3+a(1)*HRPE3)-k8*a(3)*O/(KM4+O)*a(1)/(KM5+a(1)); ...
            k4*a(1)*HRPE1/(KM3+a(1)*HRPE1)-k5*a(2); ...
            k6*a(1)*HRPE2/(KM3+a(1)*HRPE2)-k7*a(3)];
[t0,pos] = ode45(f,[0,6*3600],[0,0,0]);
RAPz = pos(end,1)
mRNAz = pos(end,2)
pco1z = pos(end,3)

%% Hypoxia simulation (O = 1)
O = 1;
f = @(t,a) [-k3*PCO4*O/(KM1+O)*a(1)/(KM2+a(1))+k1-k2*a(1)+k9*a(1)*HRPE3/(KM3+a(1)*HRPE3)-k8*a(3)*O/(KM4+O)*a(1)/(KM5+a(1)); ...
            k4*a(1)*HRPE1/(KM3+a(1)*HRPE1)-k5*a(2); ...
            k6*a(1)*HRPE2/(KM3+a(1)*HRPE2)-k7*a(3)];
[t1,posn] = ode45(f,[0,250*60],[RAPz,mRNAz,pco1z]);

%% Plot
plot(t1'/60,posn(:,2)', 'r--', 'Linewidth', 2, 'DisplayName', 'NLUC/FLUC simulation');
ylim([0 1.6])
hold on
errorbar(time/60,mRNA,sd2,'ko:' ,'Linewidth', 2, 'DisplayName', 'Data point');

plot(t1'/60,log10(posn(:,1))/2', 'b--', 'Linewidth', 2, 'DisplayName', 'RAP2.12 simulation');
ylabel('NLUC/FLUC','fontsize',14);
xlabel('Time (min)','fontsize',14);
x0=10;
y0=10;
width=550*1.4;
height=480;
legend('show','location','northwest');
legend()
set(gca,'fontsize',14)
h1=gca;
h2 = axes('Position',get(h1,'Position'));
set(h2,'YAxisLocation','right','Color','none','XTicklabel',[]);
title('UP','FontSize',14);
axis([0 120 0 1.6])
set(gca,'fontsize',14);
xticks([0 5 10 15 30 60 120 240])
