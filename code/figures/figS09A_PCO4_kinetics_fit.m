% figS09A_PCO4_kinetics_fit.m
% Reproduces Figure S9A (PCO4 kinetics fitting/visualization) from the paper.
% Robust nonlinear fit of PCO4 activity as a function of RAP2.12 and oxygen,
% then 3D surface visualization with data overlays.
% Requires Statistics and Machine Learning Toolbox (nlinfit).

clear; close all;

%% Input data (PCO4)
ERFVII = [15.625 31.25 62.5 125 250 500 1000 2000 4000 1000*ones(1,8)];
Ox     = [20*ones(1,9), 0 1.5 3 5.5 11 20 40 60];
vo     = [1.470724 2.709730 4.716129 9.452704 15.811300 24.016820 25.741070 21.844040 16.023140 ...
          1.6363 3.2727 6.5454 10.5273 19.6363 25.0909 31.6363 36.3273];

%% Robust nonlinear fitting (Welsch weights)
opts = statset('nlinfit');
opts.RobustWgtFun = 'welsch';

% Model: v = b1*(ERFVII/(b2+ERFVII))*(Ox/(b3+Ox))
% Note: keep exactly as original (ERFVII/Ox are taken from workspace)
model = @(b,x) b(1).*ERFVII./(b(2)+ERFVII).*Ox./(b(3)+Ox);

initialguess = [1 1 1];
[beta,R,J,CovB,MSE] = nlinfit([ERFVII;Ox], vo, model, initialguess, opts);

MSE
Kcat  = beta(1)/60/(1000/27311.7) % /s, MW PCO4 = 27311.7
Kmrap = beta(2)
KmO2  = beta(3)

%% Surface visualization + data overlays
set(gcf, 'Position', [100, 100, 500, 400])

ERF = linspace(0,4000);
OX  = linspace(0,60);

for i = 1:length(ERF)
    for j = 1:length(OX)
        MODEL(i,j) = beta(1)*ERF(i)/(beta(2)+ERF(i))*OX(j)/(beta(3)+OX(j));
    end
end

figure('Position', [100, 100, 800, 600]);
surf(ERF, OX, MODEL, 'EdgeColor', 'none', 'DisplayName', 'Fitting')
alpha 0.5
xlabel('RAP2.12 Concentration (µM)','FontSize',14)
ylabel('Oxygen (%)','FontSize',14)
zlabel('PCO4 activity (µmoles/min/mg)','FontSize',14)
title('PCO4','FontSize',16);
set(gca, 'FontSize', 14);
colormap hsv
colormap(gray)
shading flat
grid on
hold on

% Overlay datapoints (RAP scan at 20% O2, and O2 scan at 1250 µM RAP)
ERFVII1 = [15.625 31.25 62.5 125 250 500 1000 2000 4000];
vo1     = [1.470724 2.709730 4.716129 9.452704 15.811300 24.016820 25.741070 21.844040 16.023140];

Ox1 = [0 1.5 3 5.5 11 20 40 60];
vo2 = [1.6363 3.2727 6.5454 10.5273 19.6363 25.0909 31.6363 36.3273];

line1 = 20*ones(size(ERFVII1));
line2 = 1250*ones(size(Ox1));

plot3(ERFVII1, line1, vo1, 'r:o', 'LineWidth', 1.5, 'DisplayName', 'ERFVII data point');
plot3(line2, Ox1,  vo2, 'b:o', 'LineWidth', 1.5, 'DisplayName', 'Oxygen data point');

% Error bars (as in original)
errlRAP = [.1 .1 .1 1.3228 2.067 3.307 3.2244 1.86 2.067]/2;
errhRAP = [.1 .1 .1 1.3228 2.067 3.307 3.2244 1.86 2.067]/2;
errlOx  = [.2 .2 .2 .2 5.244 .1 4.748 3.189]/2;
errhOx  = [.2 .2 .2 .2 5.244 .1 4.748 3.189]/2;

plot3([ERFVII1',ERFVII1']', [line1',line1']', [-errlRAP',errhRAP']'+vo1, '-r','LineWidth', 1.5)
plot3([line2',line2']',    [Ox1',Ox1']',      [-errlOx',errhOx']'+vo2,    '-b','LineWidth', 1.5)
%legend('show', 'Location', 'northeast')
