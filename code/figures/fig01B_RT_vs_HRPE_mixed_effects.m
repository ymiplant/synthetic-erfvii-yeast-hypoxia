% fig01B_RT_vs_HRPE_mixed_effects.m
% Reproduces Figure 1B from the paper.
% Mixed-effects model relating response time (t90%-t10%) to HRPE density,
% controlling for promoter length (z-scored). Plots model fit with 95% CI.
% Requires Statistics and Machine Learning Toolbox (fitlme).

% ============================================
% Mixed Effects Model: Response Time ~ HRPE Density + Promoter Length (scaled)
% ============================================

clear; close all

%% Input data: response time (from Fig 1A fits) and HRPE/promoter information
geneNames = {'ADH1', 'HB1', 'HRE2', 'HUP7', 'PCO1', 'HRA1', 'LBD41', 'PDC1', 'SUS4'}';
ResponseTime = [68.54, 94.79, 57.79, 83.81, 47.62, 13.88, 14.57, 73.83, 64.31]';  % t90%-t10% (min)

% HRPE motif counts and promoter lengths (bp)
geneNames_hrpe = {'SUS4', 'HB1', 'PCO1', 'HUP7', 'HRE2', 'HRA1', 'LBD41', 'PDC1', 'ADH1'};
hrpe_number = [3, 1, 3, 5, 3, 8, 10, 4, 3];
prom_length = [2348, 518, 1212, 1614, 2402, 1052, 2262, 2632, 867];

% HRPE density (motifs per bp)
hrpe_density = hrpe_number ./ prom_length;

%% Match HRPE/promoter info to the ResponseTime gene order
hrpe_map = containers.Map(geneNames_hrpe, 1:length(geneNames_hrpe));
hrpe_density_ordered = zeros(size(ResponseTime));
prom_length_ordered = zeros(size(ResponseTime));

for i = 1:length(geneNames)
    idx = hrpe_map(geneNames{i});
    hrpe_density_ordered(i) = hrpe_density(idx);
    prom_length_ordered(i) = prom_length(idx);
end

%% Scale promoter length (z-score)
prom_length_scaled = (prom_length_ordered - mean(prom_length_ordered)) / std(prom_length_ordered);

%% Construct data table for mixed-effects model
T = table(ResponseTime, hrpe_density_ordered, prom_length_scaled, categorical(geneNames), ...
    'VariableNames', {'ResponseTime', 'HRPE_Density', 'PromLength_Scaled', 'Gene'});

%% Fit mixed-effects model (random intercept per gene)
lme = fitlme(T, 'ResponseTime ~ HRPE_Density + PromLength_Scaled + (1|Gene)');

%% Generate predictions across HRPE density (holding promoter length at mean)
xfit = linspace(min(hrpe_density_ordered)*0.9, max(hrpe_density_ordered)*1.1, 200)';

meanProm = mean(prom_length_scaled);
T_pred = table(xfit, repmat(meanProm, size(xfit)), categorical(repmat("SUS4", size(xfit))), ...
    'VariableNames', {'HRPE_Density', 'PromLength_Scaled', 'Gene'});

[yfit, yCI] = predict(lme, T_pred);

%% Compute pseudo-R^2 (training set)
fitted_train = fitted(lme);
y_true = T.ResponseTime;
Rsq = 1 - sum((y_true - fitted_train).^2) / sum((y_true - mean(y_true)).^2);

%% Plot: data, fitted relationship, and confidence interval
figure('Color', 'w', 'Position', [200, 200, 900, 650]);
hold on;

% Confidence interval band
fill([xfit; flipud(xfit)], [yCI(:,1); flipud(yCI(:,2))], ...
    [0.8 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

% Regression line
plot(xfit, yfit, 'Color', [0 0.2 0.6], 'LineWidth', 3);

% Data points
scatter(hrpe_density_ordered, ResponseTime, 120, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0 0.4470 0.7410], 'LineWidth', 1.5);

% Gene labels
for i = 1:length(geneNames)
    xOffset = 0.0003;
    yOffset = 2 * sign(ResponseTime(i) - mean(ResponseTime));
    text(hrpe_density_ordered(i) + xOffset, ResponseTime(i) + yOffset, geneNames{i}, ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'none');
end

% Axes labels
xlabel('HRPE Density (motifs / bp)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Response Time (t_{90%} - t_{10%}) [min]', 'FontSize', 14, 'FontWeight', 'bold');

% Title with model statistics (p-values taken from fixed effects coefficients)
title(sprintf(['Mixed Effects Model:\nResponse Time ~ HRPE Density + Promoter Length\n' ...
    'R^2 = %.2f, p_{HRPE} = %.1e, p_{PromLength} = %.1e'], ...
    Rsq, lme.Coefficients.pValue(2), lme.Coefficients.pValue(3)), ...
    'FontSize', 16, 'FontWeight', 'bold');

% Style
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
box on;
grid on;
xlim([min(hrpe_density_ordered)*0.85, max(hrpe_density_ordered)*1.15]);
ylim([min(ResponseTime)*0.85, max(ResponseTime)*1.15]);
hold off;
