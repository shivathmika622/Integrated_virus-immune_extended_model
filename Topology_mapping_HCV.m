%% 1. Load the Data
load('Feasible_Data_For_Clustering.mat', 'dataNormalized', 'feasibleData');

%% 2. Run SVD (Singular Value Decomposition)
% SVD is the mathematical engine behind PCA.
% U: Left singular vectors (related to scores)
% S: Singular values (related to variance)
% V: Right singular vectors (related to coefficients/loadings)
[U, S, V] = svd(dataNormalized, 'econ');

%% 3. Calculate Scores and Variance
% Scores = U * S
score = U * S; 

% Variance explained calculation
latent = diag(S).^2 / (size(dataNormalized,1) - 1);
varExplained = 100 * latent / sum(latent);
fprintf('First two components explain %.2f%% of the variance.\n', sum(varExplained(1:2)));

%% 4. Plot the PCA Projection
figure;
% score(:,1) is PC1, score(:,2) is PC2
scatter(score(:,1), score(:,2), 60, log10(feasibleData.VT_96h), 'filled', 'MarkerFaceAlpha', 0.7);
grid on;
h = colorbar;
ylabel(h, 'Log10(Steady-State Titer V_T)');
colormap(jet);

xlabel(['PC1 (', num2str(varExplained(1), '%.1f'), '% Variance)']);
ylabel(['PC2 (', num2str(varExplained(2), '%.1f'), '% Variance)']);
title('Virus Fitness Manifold (SVD-based PCA)');

%% 5. Identify Driving Parameters
% V matrix carries the weights of each parameter
paramNames = feasibleData.Properties.VariableNames(1:20);
[~, idx] = sort(abs(V(:,1)), 'descend');

fprintf('\nTop 3 Parameters driving the separation (PC1):\n');
for j = 1:3
    fprintf('%d. %s (Weight: %.3f)\n', j, paramNames{idx(j)}, V(idx(j),1));
end