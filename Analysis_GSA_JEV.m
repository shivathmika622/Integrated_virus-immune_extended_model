clear; close all; clc;

%% 1. Load Results
Results = readtable('GSA_JEV_Results.csv');

all_data = Results(~isnan(Results.VT_96h), :);

figure;
histogram(log10(all_data.VT_96h + 1e-6),50);
xlabel('log_{10}(VT_{96h})');
ylabel('Count');
title('Distribution of Viral Titer at 96 h');

%% 2. All parameter names
paramNames = { ...
    'k_en', 'k_a', 'k_r', 'k_t', 'k_c', 'rcsat', 'tau', ...
    'k_RIGI', 'k_MAVS', 'k72', 'mu_ISG_RNA', 'mu_p', ...
    'M', 'Omega', 'gamma_RIGI', 'threshold', 'n' ...
};

%% 3. PRCC Analysis
clean_data = Results(~isnan(Results.VT_96h), :);
X_all = table2array(clean_data(:, paramNames));
Y_output = clean_data.VT_96h;

data_combined = [X_all, Y_output];
[N_samples, N_vars] = size(data_combined);
Ranks = zeros(N_samples, N_vars);
for j = 1:N_vars
    [~, ~, Ranks(:,j)] = unique(data_combined(:,j));
end
C = corrcoef(Ranks);
C_inv = inv(C);

rho = zeros(1, N_vars-1);
p_strings = cell(N_vars-1, 1);
df = N_samples - N_vars;

for j = 1:N_vars-1
    rho(j) = -C_inv(j, N_vars) / sqrt(C_inv(j,j) * C_inv(N_vars, N_vars));
    t_stat = rho(j) * sqrt(df / (1 - rho(j)^2));
    p_val = 2 * (1 - tcdf(abs(t_stat), df));
    if p_val < 1e-16
        p_strings{j} = '< 1.0e-16';
    elseif p_val < 0.001
        p_strings{j} = '< 0.001';
    else
        p_strings{j} = sprintf('%.5f', p_val);
    end
end

PRCC_Results = table(paramNames', rho', p_strings, ...
    'VariableNames', {'Parameter', 'PRCC_Coefficient', 'P_Value'});
disp('--- PRCC Results ---'); disp(PRCC_Results);

%% 4. PRCC Plot
figure('Position', [100, 100, 650, 400]);
bar(rho, 'FaceColor', [0.2 0.6 0.4]);
set(gca, 'XTick', 1:length(paramNames), 'XTickLabel', paramNames, ...
    'TickLabelInterpreter', 'none', 'FontSize', 9);
xtickangle(45);
ylabel('PRCC Coefficient', 'FontWeight', 'bold');
title('GSA: PRCC Impact on Viral Titer (96h)');
grid on;

%% 5. Filter Strong Parameters (|PRCC| > 0.15)
strong_idx = abs(rho) > 0.15;
strongParams = paramNames(strong_idx);
fprintf('\nStrong parameters selected: %d out of %d\n', sum(strong_idx), length(paramNames));
disp(strongParams');

%% 6. Feasible Data
all_data = Results(~isnan(Results.VT_96h), :);
feasible_data = all_data(all_data.Feasible == 1, :);
fprintf('Feasible samples: %d (%.1f%%)\n', ...
    height(feasible_data), ...
    100*height(feasible_data)/height(all_data));
writetable(feasible_data, 'GSA_JEV_Only_Feasible_Profiles.csv');

%% 7. Prepare data — strong params only
% Use linear data to maintain correct variance; DO NOT take log10 of inputs.
X_allparams = table2array(all_data(:, paramNames));
X_scaled = zscore(X_allparams);
fprintf('Data matrix size for PCA/tSNE: %d x %d\n', size(X_scaled));

%% 8. PCA — on all params
fprintf('Computing PCA...\n');
[coeff, ~, ~, ~, explained] = pca(X_scaled);

fprintf('PC1: %.1f%% | PC2: %.1f%% | Total(PC1+PC2): %.1f%%\n', ...
    explained(1), explained(2), explained(1)+explained(2));

%% 9. t-SNE
fprintf('Computing t-SNE...\n');
% Swapped 'exact' to 'barneshut' to completely eliminate Out of Memory errors.
X_tsne = tsne(X_scaled, 'Algorithm', 'barneshut', 'Distance', 'euclidean', ...
    'NumDimensions', 2, 'Perplexity', 150, 'Exaggeration', 12, ...
    'Standardize', false);

%% 10. K-Means (K=2)
fprintf('K-Means clustering...\n');
opts_km = statset('Display', 'final', 'MaxIter', 500);
[cluster_assignments, ~] = kmeans(X_scaled, 2, ...
    'Distance', 'sqeuclidean', 'Replicates', 20, 'Options', opts_km);

%% 11. Cluster centroids back-projection
n_params = length(paramNames);
centroids = zeros(2, n_strong);
for k = 1:2
    centroids(k,:) = mean(X_strong(cluster_assignments == k, :), 1);
end

Strategy_Profiles = table(ParamNames', ...
    (centroids(1,:))', ...
    (centroids(2,:))', ...
    'VariableNames', {'Parameter', 'Cluster1_Mean', 'Cluster2_Mean'});
disp('--- Strategy Mapping Table ---'); disp(Strategy_Profiles);

%% 12. Plots
figure('Position', [50, 150, 1150, 380]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');

% Left: t-SNE
nexttile;
gscatter(X_tsne(:,1), X_tsne(:,2), cluster_assignments, ...
    [0.9 0.0 0.9; 0.0 0.9 0.9], '.', 6);
set(gca, 'FontSize', 9);
xlabel('t-SNE Dim 1'); ylabel('t-SNE Dim 2');
title('Viral Manifold Topology', 'FontWeight', 'bold');
legend('Cluster 1', 'Cluster 2', 'Location', 'best', 'FontSize', 8);
grid on;

% Middle: PC1 loadings
nexttile;
bar(coeff(1:n_strong, 1), 'FaceColor', [0.1 0.5 0.7]);
set(gca, 'XTick', 1:n_strong, 'XTickLabel',ParamNames, ...
    'TickLabelInterpreter', 'none', 'FontSize', 8.5);
xtickangle(60);
ylabel('PC1 Loading');
title(['PC1 (', num2str(explained(1),'%.1f'), '% Var)'], 'FontWeight', 'bold');
grid on;

% Right: PC2 loadings
nexttile;
bar(coeff(1:n_strong, 2), 'FaceColor', [0.6 0.2 0.5]);
set(gca, 'XTick', 1:n_strong, 'XTickLabel',ParamNames, ...
    'TickLabelInterpreter', 'none', 'FontSize', 8.5);
xtickangle(60);
ylabel('PC2 Loading');
title(['PC2 (', num2str(explained(2),'%.1f'), '% Var)'], 'FontWeight', 'bold');
grid on;
