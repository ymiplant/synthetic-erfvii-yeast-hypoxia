% fig06_GG.m
% Reproduces Figure 6 (GG model) from the paper.
% Note: this script also generates the RAP2.12 fitting plot (related to Fig S9B),
% then generates the GG NLUC/FLUC simulation overlapped with data.

clear; close all;

%% Part 1: RAP2.12 fitting plot (related to Fig S9B)
t = [0 0.5*3600 1*3600 3*3600];
RAP = 150*[0.1 0.188053783093125 0.300356526688525 0.938186923750606]; %add t0=0
sd1=150*[0.001 0.022762065 0.047373145 0.112143258];
k1 = 0.0301;
k2 = 0.0001; 
k3 = 26.5878; %per second
KM1 = 16.4217; %percent
KM2 = 270.2505; %micro-molar
O = 21;
RAP0 = RAP(1);
PCO4 = 0.03;

O = 1;

f = @(t,a) [-k3*PCO4*O/(KM1+O)*a(1)/(KM2+a(1))+k1-k2*a(1)];
[t1,posn] = ode45(f,[0,250*60],[RAP0]);
plot(t1'/60,posn(:,1)', 'k--', 'Linewidth', 2, 'DisplayName', 'Simulation')
ylabel('Adjusted RAP2.12 (µM)','fontsize',14);
xlabel('Time (min)','fontsize',14);
hold on
errorbar(t/60,RAP,sd1,'ko-' ,'Linewidth', 2, 'DisplayName', 'Data point')
legend('show','location','northwest');
legend()

set(gca,'fontsize',14)
h1=gca;
h2 = axes('Position',get(h1,'Position'));
set(h2,'YAxisLocation','right','Color','none','XTicklabel',[]);
ylabel('NLUC/FLUC','fontsize',14);

x0=10;
y0=10;
width=550*1.5;
height=480;

set(gca,'fontsize',14);

%% Part 2: GG model NLUC/FLUC simulation (Figure 6)
figure
time = [0 5*60 10*60 15*60 30*60 60*60 120*60 240*60];
mRNA = [0.229763793	0.330703629	0.339784906	0.339279166	0.323545532	0.338143938	0.359982939	0.363713039];
sd2=[0.006966078 0.04362041 0.037627681	0.049687797	0.063665912	0.052152682	0.02492597 0.053209287];

mRNA0 = mRNA(1);
O = 21;
alpha = -k3*PCO4*O/(KM1+O);
RAPz = roots([-k2,-k2*KM2+k1+alpha,k1*KM2]);
RAP0 = RAPz(RAPz>0);
HRPE1 = 0.07; %
k4 = 0.006;
k5 = 0.014;
KM3 = 0.5;
PCO4 = 0.06;

f = @(t,a) [-k3*PCO4*O/(KM1+O)*a(1)/(KM2+a(1))+k1-k2*a(1); ...
            k4*a(1)*HRPE1/(KM3+a(1)*HRPE1)-k5*a(2)];
[t0,pos] = ode45(f,[0,6*3600],[0.1,0]);
RAPz = pos(end,1)
mRNAz = pos(end,2)

O = 1;
f = @(t,a) [-k3*PCO4*O/(KM1+O)*a(1)/(KM2+a(1))+k1-k2*a(1); ...
            k4*a(1)*HRPE1/(KM3+a(1)*HRPE1)-k5*a(2)];
[t1,posn] = ode45(f,[0,250*60],[RAPz,mRNAz]);

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
title('GG','FontSize',14);
set(gca,'fontsize',14)
h1=gca;
h2 = axes('Position',get(h1,'Position'));
set(h2,'YAxisLocation','right','Color','none','XTicklabel',[]);
ylim([0 1.6])

set(gca,'fontsize',14);
