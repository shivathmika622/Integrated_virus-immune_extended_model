%% 1. Load the Results
resultsTable = readtable('GSA_HCV_Results.csv'); 

%% 2. Filter for Feasible Sets
feasibleData = resultsTable(resultsTable.Feasible == 1, :);

% Counts display
fprintf('Total Simulations: %d\n', height(resultsTable));
fprintf('Feasible Strategies Found: %d\n', height(feasibleData));

%% 3. Manual Normalization (Alternative to zscore)
% Parameter columns (1 to 20)
paramColumns = resultsTable.Properties.VariableNames(1:20);
dataMatrix = table2array(feasibleData(:, paramColumns));

% Log transformation
logData = log10(dataMatrix);

% Manual Z-score calculation: (x - mean) / std
mu = mean(logData);
sigma = std(logData);
dataNormalized = (logData - mu) ./ sigma;

%% 4. Save Feasible Cloud
writetable(feasibleData, 'Feasible_Manifold_Points.csv');
save('Feasible_Data_For_Clustering.mat', 'dataNormalized', 'feasibleData');

disp('Step 1 Completed: Feasible cloud is ready for Visualization.');