% figS09A_PCO1_kinetics_fit.m
% Reproduces Figure S9A (PCO1 kinetics fitting/visualization) from the paper.
% Robust nonlinear fit of PCO1 activity as a function of RAP2.12 and oxygen,
% then 3D surface visualization with data overlays.
% Requires Statistics and Machine Learning Toolbox (nlinfit).

clear; close all;

%% Input data (PCO1)
ERFVII = [15.625 31.25 62.5 125 250 500 1000 2000 4000 1250*ones(1,8)];
Ox     = [20*ones(1,9) 0 1.5 3 5.5 11 20 40 60];
vo     = [.339 .756 1.549 2.264 3.527 5.205 5.859 6.244 6.657 .5454 1.6545 2.5818 4.1273 6.1091 9.3273 11.0909 13.2364];

%% Robust nonlinear fitting (Talwar weights)
opts = statset('nlinfit');
opts.RobustWgtFun = 'Talwar';

% Model: v = b1*(ERFVII/(b2+ERFVII))*(Ox/(b3+Ox))
% Note: keep exactly as original (ERFVII/Ox are taken from workspace)
model = @(b,x) b(1).*ERFVII./(b(2)+ERFVII).*Ox./(b(3)+Ox);

initialguess = [1 1 1];
[beta,R,J,CovB,MSE] = nlinfit([ERFVII;Ox], vo, model, initialguess, opts);

MSE
Kcat  = beta(1)/60/(1000/33013.4) % /s, MW PCO1 = 33013.4
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
zlabel('PCO1 activity (µmoles/min/mg)','FontSize',14)
title('PCO1','FontSize',16);
set(gca, 'FontSize', 14);
colormap hsv
colormap(gray)
shading flat
grid on
hold on

% Overlay datapoints (RAP scan at 20% O2, and O2 scan at 1250 µM RAP)
ERFVII1 = [15.625 31.25 62.5 125 250 500 1000 2000 4000];
vo1     = [.339 .756 1.549 2.264 3.527 5.205 5.859 6.244 6.657];

Ox1 = [0 1.5 3 5.5 11 20 40 60];
vo2 = [.5454 1.6545 2.5818 4.1273 6.1091 9.3273 11.0909 13.2364];

line1 = 20*ones(size(ERFVII1));
line2 = 1250*ones(size(Ox1));

plot3(ERFVII1, line1, vo1, 'r:o', 'LineWidth', 1.5, 'DisplayName', 'ERFVII data point');
plot3(line2, Ox1,  vo2, 'b:o', 'LineWidth', 1.5, 'DisplayName', 'Oxygen data point');

% Error bars (as in original)
errlRAP = [.1 .1 .1 .1 .3 .59 .7 .59 1.016]/2;
errhRAP = [.1 .1 .1 .1 .3 .59 .7 .59 1.016]/2;
errlOx  = [.1 .1 .1 .1 .1 .1 .886 2.95]/2;
errhOx  = [.1 .1 .1 .1 .1 .1 .886 2.95]/2;

plot3([ERFVII1',ERFVII1']', [line1',line1']', [-errlRAP',errhRAP']'+vo1, '-r','LineWidth', 1.5)
plot3([line2',line2']',    [Ox1',Ox1']',      [-errlOx',errhOx']'+vo2,    '-b','LineWidth', 1.5)
