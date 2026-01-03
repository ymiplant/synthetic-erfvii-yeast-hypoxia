% fig01A_logistic_fit.m
% Reproduces Figure 1A from the paper.
% Logistic fitting of hypoxia-responsive gene expression time courses.
% Response time is defined as t90% - t10% = 2*log(9)/k.
% Requires MATLAB Optimization Toolbox (lsqcurvefit).

%% Input data
% Gene expression means and standard deviations at each time point
genes = {
    'ADH1',  [0.78, 0.97, 1.54, 4.46, 22.37, 90.27, 309.23], [0.15, 0.37, 0.34, 0.84, 4.12, 8.39, 39.29];
    'HB1',   [0.67, 1.56, 3.41, 5.12, 12.05, 27.58, 74.33], [0.26, 0.28, 0.73, 0.92, 3.04, 2.40, 17.91];
    'HRE2',  [0.86, 1.63, 4.08, 10.37, 42.45, 142.28, 218.53], [0.24, 0.53, 0.78, 1.75, 7.15, 12.07, 29.42];
    'HUP7',  [0.53, 0.93, 1.85, 3.23, 7.41, 26.13, 116.61], [0.34, 0.18, 0.58, 0.71, 2.03, 4.78, 18.16];
    'PCO1',  [0.81, 2.91, 9.68, 24.89, 64.66, 105.57, 136.63], [0.17, 0.29, 3.84, 3.77, 9.25, 11.62, 12.67];
    'HRA1',  [0.92, 15.40, 54.54, 119.91, 158.77, 171.26, 183.97], [0.08, 3.95, 7.14, 8.67, 13.85, 30.52, 21.75];
    'LBD41', [1.07, 33.24, 92.24, 173.83, 218.37, 218.86, 249.90], [0.22, 6.57, 13.54, 14.48, 35.70, 20.79, 26.40];
    'PDC1',  [0.76, 1.34, 3.89, 9.86, 39.59, 81.60, 143.80], [0.20, 0.61, 0.79, 1.64, 3.95, 13.92, 21.83];
    'SUS4',  [0.60, 0.71, 1.29, 2.54, 11.37, 45.06, 101.60], [0.27, 0.06, 0.26, 0.58, 1.33, 7.16, 27.63];
};

%% Time points (minutes after shift to hypoxia)
time = [0, 5, 10, 15, 30, 60, 120];

%% Preallocation for summary values
geneNames = cell(size(genes,1),1);
t90_t10_values = zeros(size(genes,1),1);

%% Loop over genes and fit logistic model
for i = 1:size(genes, 1)

    % Extract data for current gene
    geneName = genes{i, 1};
    expr     = genes{i, 2};
    stdev    = genes{i, 3};

    % Logistic model: A / (1 + exp(-k*(t - t_half)))
    logistic = @(params, t) params(1) ./ ...
        (1 + exp(-params(2) * (t - params(3))));

    % Initial parameter guess
    initialGuess = [max(expr), 0.05, 30];
    if strcmp(geneName, 'HB1')
        initialGuess(2) = 0.1;
    end

    % Nonlinear least-squares fit
    fitParams = lsqcurvefit(logistic, initialGuess, time, expr);

    % Extract fitted parameters
    A      = fitParams(1);
    k      = fitParams(2);
    t_half = fitParams(3);

    % Compute response time (t90% - t10%)
    t90_t10 = 2 * log(9) / k;

    % Evaluate fitted curve for plotting
    t_fit    = linspace(min(time), max(time), 200);
    expr_fit = logistic(fitParams, t_fit);

    %% Plot experimental data and fitted curve
    figure;
    errorbar(time, expr, stdev, 'ko', ...
        'LineWidth', 1.5, 'MarkerSize', 6); hold on;
    plot(t_fit, expr_fit, 'r-', 'LineWidth', 2);
    xlabel('Time (min)');
    ylabel('Gene Expression');
    title(['Gene Induction - ', geneName]);
    legend('Data ± SD', 'Logistic Fit', 'Location', 'NorthWest');
    grid on;

    % Annotate response time on the plot
    xlim_vals = xlim;
    ylim_vals = ylim;
    yRange = ylim_vals(2) - ylim_vals(1);
    maxExpr = max([expr, expr_fit]);
    text_y = maxExpr + 0.05 * yRange;

    if text_y > ylim_vals(2)
        ylim([ylim_vals(1), text_y + 0.05 * yRange]);
    end

    text_x = xlim_vals(1) + 0.65 * (xlim_vals(2) - xlim_vals(1));
    text(text_x, text_y, ...
        ['$$t_{90\%} - t_{10\%} = ', ...
         num2str(t90_t10, '%.2f'), '\ \mathrm{min}$$'], ...
        'Interpreter', 'latex', 'FontSize', 12, ...
        'FontWeight', 'bold', 'Color', 'blue');

    %% Store results and print summary
    geneNames{i} = geneName;
    t90_t10_values(i) = t90_t10;

    fprintf('\n[%s] Fitted Parameters:\n', geneName);
    fprintf('  Max expression (A)         = %.3f\n', A);
    fprintf('  Rate constant (k)          = %.4f\n', k);
    fprintf('  Half-max time (t_1/2)      = %.2f min\n', t_half);
    fprintf('  t90%% - t10%% response time = %.2f min\n', t90_t10);
end

%% Summary table printed to console
fprintf('\nSummary of t90%%-t10%% response time for all genes:\n');
fprintf('Gene\tResponse time (min)\n');
for i = 1:length(geneNames)
    fprintf('%s\t%.2f\n', geneNames{i}, t90_t10_values(i));
end

