% fig06_UG.m
% Reproduces Figure 6 (UG model) from the paper.
% ODE simulation of NLUC/FLUC output for the UG circuit, overlapped with data.

clear; close all;

%% Experimental data (NLUC/FLUC)
time = [0 5*60 10*60 15*60 30*60 60*60 120*60 4*3600];
mRNA = [0.225138928	0.318223499	0.616141415	0.599652336	0.600744775	0.712851459	0.585360651	0.576052281]; %
sd2=[0.029499742 0.155438362 0.123617956 0.258974589 0.111630525 0.06084104	0.035441158	0.018654601];

%% Parameters (as used in the paper)
k1 = 0.0301;
k2 = 0.0001;
k3 = 26.5878; %per second
KM1 = 16.4217; %percent
KM2 = 270.2505; %micro-molar
k4 = 0.0067;
k5 = 0.011;
KM3 = 0.5;
PCO4 = 0.08;
k9 = 0.05;
HRPE1 = 0.04; %
HRPE3 = 0.02; %

%% Initial conditions under normoxia (O = 21)
mRNA0 = mRNA(1);
O = 21;

RAP0 = k5*mRNA0*KM3/(k4*HRPE1-HRPE1*k5*mRNA0);
X = k3*PCO4*O/(KM1+O)*RAP0/(KM2+RAP0)-k1+k2*RAP0;

f = @(t,a) [-k3*PCO4*O/(KM1+O)*a(1)/(KM2+a(1))+k1-k2*a(1)+k9*a(1)*HRPE3/(KM3+a(1)*HRPE3); ...
            k4*a(1)*HRPE1/(KM3+a(1)*HRPE1)-k5*a(2)];
[t0,pos] = ode45(f,[0,6*3600],[0.1,0]);
RAPz = pos(end,1)
mRNAz = pos(end,2)

%% Hypoxia simulation (O = 1)
O = 1;
f = @(t,a) [-k3*PCO4*O/(KM1+O)*a(1)/(KM2+a(1))+k1-k2*a(1)+k9*a(1)*HRPE3/(KM3+a(1)*HRPE3); ...
            k4*a(1)*HRPE1/(KM3+a(1)*HRPE1)-k5*a(2)];
[t1,posn] = ode45(f,[0,250*60],[RAPz,mRNAz]);

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
legend('show','location','northwest', 'box' , 'on' );
legend()
title('UG','FontSize',14);
set(gca,'fontsize',14)
h1=gca;
h2 = axes('Position',get(h1,'Position'));
set(h2,'YAxisLocation','right','Color','none','XTicklabel',[]);
ylim([0 1.6])
set(gca,'fontsize',14);
