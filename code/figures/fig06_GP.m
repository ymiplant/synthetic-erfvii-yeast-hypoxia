% fig06_GP.m
% Reproduces Figure 6 (GP model) from the paper.
% ODE simulation of NLUC/FLUC output for the GP circuit, overlapped with data.

clear; close all;

%% Experimental data (NLUC/FLUC)
time=[0 5*60 10*60 15*60 30*60 60*60 120*60 240*60];
mRNA=[0.154618429 0.347964663 0.399412914 0.393142508 0.371229577 0.449077333 0.300324936 0.369111386];
sd2=[0.049263415 0.153708535 0.040567688 0.052965391 0.019806177 0.123992481 0.069462394 0.049302054];

%% Parameters (as used in the paper)
k1 = 0.0301;
k2 = 0.0001; 
k3 = 26.5878; %Kcat for PCO4 per second
KM1 = 16.4217; %KM O2 for PCO4 percent
KM2 = 270.2505; %KM rap for PCO4 micro-molar
k4 = 0.0067;
k5 = 0.012;
KM3 = 0.5;
PCO4 = 0.06;
k6 = 0.003;
k7 = 0.0009;
k8 = 9.6688; % Kcat for PCO1
KM4 = 19.5138; % KM O2 for PCO1
KM5 = 286.6942; % KM rap for PCO1
HRPE1 = 0.06; %
HRPE2 = 0.01;

%% Initial conditions under normoxia (O = 21)
mRNA0 = mRNA(1);
O = 21;

RAP0 = k5*mRNA0*KM3/(k4*HRPE1-HRPE1*k5*mRNA0);

f = @(t,a) [-k3*PCO4*O/(KM1+O)*a(1)/(KM2+a(1))+k1-k2*a(1)-k8*a(3)*O/(KM4+O)*a(1)/(KM5+a(1)); ...
            k4*a(1)*HRPE1/(KM3+a(1)*HRPE1)-k5*a(2); ...
            k6*a(1)*HRPE2/(KM3+a(1)*HRPE2)-k7*a(3)];
[t0,pos] = ode45(f,[0,6*3600],[0.1,0,0]);
RAPz = pos(end,1)
mRNAz = pos(end,2)
pco1z = pos(end,3)

%% Hypoxia simulation (O = 1)
O = 1;
f = @(t,a) [-k3*PCO4*O/(KM1+O)*a(1)/(KM2+a(1))+k1-k2*a(1)-k8*a(3)*O/(KM4+O)*a(1)/(KM5+a(1)); ...
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
width=550*1.5;
height=480;
legend('show','location','northwest');
legend()
title('GP','FontSize',14);
set(gca,'fontsize',14)
h1=gca;
h2 = axes('Position',get(h1,'Position'));
set(h2,'YAxisLocation','right','Color','none','XTicklabel',[]);
ylim([0 1.6])
set(gca,'fontsize',14);
