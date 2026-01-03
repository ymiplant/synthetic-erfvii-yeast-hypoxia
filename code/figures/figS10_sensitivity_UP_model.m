% figS10_sensitivity_UP_model.m
% Reproduces Figure S10 (UP model sensitivity analysis) from the paper.
% Single-parameter sensitivity scan (80–120% of fitted value) for UP model parameters.
% Outputs:
%   - Multi-panel curves: t90%-t10% vs parameter scaling (S10C-style)
%   - Ranking bar plot: parameters sorted by response-time sensitivity (S10B-style)

clear; close all; clc;

%% Base parameters for UP model
base = struct(...
    'k1', 0.0301, ...
    'k2', 0.0001, ...
    'k3', 26.5878, ...
    'KM1', 16.4217, ...
    'KM2', 270.2505, ...
    'PCO4', 0.08, ...
    'k4', 0.0088, ...
    'k5', 0.012, ...
    'k6', 0.001, ...
    'k7', 0.0009, ...
    'k8', 9.6688, ...
    'KM3', 0.5, ...
    'KM4', 19.5138, ...
    'KM5', 286.6942, ...
    'k9', 0.04, ...
    'HRPE1', 0.02, ...
    'HRPE2', 0.005, ...
    'HRPE3', 0.015, ...
    'O2', 1);

%% Parameters to scan (80–120% of base)
param_names = {'k1','k2','k3','k4','k5','k6','k7','k8','k9', ...
               'KM1','KM2','KM3','KM4','KM5', ...
               'HRPE1','HRPE2','HRPE3'};
n_param = numel(param_names);

scan_list = cell(n_param, 2);
for i = 1:n_param
    pname = param_names{i};
    base_val = base.(pname);
    scan_list{i,1} = pname;
    scan_list{i,2} = linspace(0.8*base_val, 1.2*base_val, 20);
end

%% Sensitivity computation: response time (t90%-t10%) for each parameter value
t90_t10_all = zeros(n_param, 1);
t90_t10_curves = cell(n_param, 1);
global_min = inf;
global_max = -inf;

for p = 1:n_param
    pname = scan_list{p,1};
    pvals = scan_list{p,2};
    tdiff = zeros(size(pvals));

    for i = 1:length(pvals)
        param = base;
        param.(pname) = pvals(i);

        % --- Normoxia pre-equilibration (O2 = 20%) to get initial condition
        param.O2 = 20;
        f_norm = @(t,a) [
            param.k1 ...
            - param.k2*a(1) ...
            - param.k3*param.PCO4*param.O2/(param.KM1+param.O2)*a(1)/(param.KM2+a(1)) ...
            - param.k8*a(3)*param.O2/(param.KM4+param.O2)*a(1)/(param.KM5+a(1)) ...
            + param.k9*a(1)*param.HRPE1/(param.KM3 + a(1)*param.HRPE1);

            param.k4*a(1)*param.HRPE3/(param.KM3 + a(1)*param.HRPE3) - param.k5*a(2);

            param.k6*a(1)*param.HRPE2/(param.KM3 + a(1)*param.HRPE2) - param.k7*a(3)
        ];

        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [~, y_norm] = ode45(f_norm, [0, 6*3600], [0,0,0], options);
        init = y_norm(end,:);

        % --- Hypoxia simulation (O2 = 1%)
        param.O2 = 1;
        f_hypo = @(t,a) [
            param.k1 ...
            - param.k2*a(1) ...
            - param.k3*param.PCO4*param.O2/(param.KM1+param.O2)*a(1)/(param.KM2+a(1)) ...
            - param.k8*a(3)*param.O2/(param.KM4+param.O2)*a(1)/(param.KM5+a(1)) ...
            + param.k9*a(1)*param.HRPE1/(param.KM3 + a(1)*param.HRPE1);

            param.k4*a(1)*param.HRPE3/(param.KM3 + a(1)*param.HRPE3) - param.k5*a(2);

            param.k6*a(1)*param.HRPE2/(param.KM3 + a(1)*param.HRPE2) - param.k7*a(3)
        ];

        [t_hypo, y_hypo] = ode45(f_hypo, [0, 10*3600], init, options);
        y_mRNA = y_hypo(:,2);

        % Normalize to [0,1] using initial value and steady-state (mean last 10 points)
        ss_val_max = mean(y_mRNA(end-9:end));
        ss_val_min = mean(y_mRNA(1));
        y_norm_mRNA = (y_mRNA - ss_val_min) / (ss_val_max - ss_val_min);

        % t10 and t90 crossing (first time reaching 10% and 90%)
        i10 = find(y_norm_mRNA >= 0.1, 1, 'first');
        i90 = find(y_norm_mRNA >= 0.9, 1, 'first');

        if isempty(i10) || isempty(i90)
            tdiff(i) = NaN;
        else
            t10 = t_hypo(i10);
            t90 = t_hypo(i90);
            tdiff(i) = (t90 - t10) / 60; % minutes
        end
    end

    % Sensitivity metric used for ranking (range across scan)
    t90_t10_all(p) = max(tdiff) - min(tdiff);
    t90_t10_curves{p} = tdiff;

    global_min = min(global_min, min(tdiff));
    global_max = max(global_max, max(tdiff));
end

%% Multi-panel sensitivity curves (S10C-style)
figure('Position', [100, 100, 1600, 1000]);
for p = 1:n_param
    pname = scan_list{p,1};
    pvals = scan_list{p,2};
    tdiff = t90_t10_curves{p};
    x_rel = pvals / base.(pname);

    subplot(5, 4, p); hold on;
    x_fit = linspace(min(x_rel), max(x_rel), 200);
    y_fit = spline(x_rel, tdiff, x_fit);
    plot(x_fit, y_fit, 'b-', 'LineWidth', 2);
    plot(x_rel, tdiff, 'bo', 'MarkerSize', 5, 'LineWidth', 1.5);

    title([pname, '  \Delta = ', num2str(t90_t10_all(p), '%.2f')], 'FontSize', 11);
    xlabel([pname, ' (% of base)'], 'FontSize', 10);
    ylabel('t_{90%} - t_{10%} (min)', 'FontSize', 10);
    ylim([global_min - 1, global_max + 1]);
    xlim([0.8 1.2]);
    xticks([0.8 0.9 1.0 1.1 1.2]);
    xticklabels({'80%', '90%', '100%', '110%', '120%'});
    grid on;
end
sgtitle('Sensitivity of Response Speed (t_{90%} - t_{10%}) to Model Parameters', 'FontSize', 15);

%% Ranking bar plot (S10B-style)
figure;
[sorted_delta, idx] = sort(t90_t10_all, 'descend');
sorted_names = param_names(idx);
barh(flip(sorted_delta), 'FaceColor', [0.3, 0.5, 0.85]);
set(gca, 'YTick', 1:n_param, 'YTickLabel', flip(sorted_names), 'FontSize', 12);
xlabel('\Delta Response Time (t_{90%} - t_{10%}) (min)');
title('Sensitivity of Response Speed to Parameters');
grid on;
