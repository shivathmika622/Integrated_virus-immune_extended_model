%% 1. Load Data
load('Feasible_Data_For_Clustering.mat', 'dataNormalized', 'feasibleData');
[U, S, V] = svd(dataNormalized, 'econ');
score = U * S;
X = dataNormalized; 

%% 2. Manual K-means Implementation
k = 2; 
max_iter = 100;
[N, D] = size(X);
randIdx = randperm(N, k);
centroids = X(randIdx, :);
idx = zeros(N, 1);

for iter = 1:max_iter
    old_idx = idx;
    for i = 1:N
        dist = sum((X(i,:) - centroids).^2, 2);
        [~, idx(i)] = min(dist);
    end
    if isequal(idx, old_idx), break; end
    for j = 1:k
        centroids(j,:) = mean(X(idx == j, :), 1);
    end
end

%% 3. Basic Scatter Visualization (No gscatter)
figure; hold on;
colors = ['r', 'b']; % Red for Cluster 1, Blue for Cluster 2
markers = ['o', 'x'];
for j = 1:k
    cluster_points = score(idx == j, :);
    scatter(cluster_points(:,1), cluster_points(:,2), 60, colors(j), markers(j), 'LineWidth', 1.5);
end
grid on;
xlabel('PC1'); ylabel('PC2');
legend('Cluster 1', 'Cluster 2');
title('Virus Strategy Clusters (Manual K-means)');
hold off;

%% 4. Strategy Analysis
paramNames = feasibleData.Properties.VariableNames(1:20);
fprintf('\n--- Cluster Analysis ---\n');
for j = 1:k
    cluster_pts = feasibleData(idx == j, :);
    fprintf('\nCluster %d (Size: %d):\n', j, height(cluster_pts));
    fprintf('  Avg V_T: %.2f\n', mean(cluster_pts.VT_96h));

    % Find parameter with highest variation from global mean
    cluster_mean_log = mean(log10(table2array(cluster_pts(:, 1:20))));
    global_mean_log = mean(log10(table2array(feasibleData(:, 1:20))));
    [~, top_p] = max(abs(cluster_mean_log - global_mean_log));

    fprintf('  Key Strategy Driver: %s\n', paramNames{top_p});
end