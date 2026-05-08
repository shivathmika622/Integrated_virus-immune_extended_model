clear; close all; clc;

%% Load GSA results
Results = readtable('GSA_JEV_Results.csv');

%% Parameters to analyze (20 names)
paramNames = {
    'ktV', 'kcV', 'krV', 'NCV', 'tau', ...
    'kenV', 'kaV', 'kMAVS', 'kinhISGn', 'muISGRNA', ...
    'ktISGRNA', 'M', 'Omega', 'gammaRIGI', 'kRIGI', ...
    'k72', 'degARCISGn', 'mu_pV', 'threshold', 'nHill'
    };

%% Correlation-based GSA
Sensitivity = zeros(length(paramNames),1);
for i = 1:length(paramNames)
    x = Results.(paramNames{i});
    y = Results.VT_96h;
    R = corrcoef(x,y,'Rows','complete');
    Sensitivity(i) = R(1,2);
end

% Create Table with Transpose (') to fix row error
GSATable = table(paramNames', Sensitivity, ...
    'VariableNames', {'Parameter','Sensitivity'});
disp('--- Simple Correlation Results ---')
disp(GSATable)

%% PRCC Analysis Phase
%% PRCC Analysis Phase (PURE MANUAL - No Toolbox Required)
% 1. Remove rows where ODE failed
clean_data = Results(~isnan(Results.VT_96h), :); 

% 2. Extract Data
X = table2array(clean_data(:, paramNames)); 
Y_output = clean_data.VT_96h;

% 3. Manual Rank Transformation (Standard logic)
data_combined = [X, Y_output];
[N_samples, N_vars] = size(data_combined);
Ranks = zeros(N_samples, N_vars);

for j = 1:N_vars
    [~, ~, Ranks(:,j)] = unique(data_combined(:,j));
end

% 4. Manual Correlation Matrix using 'corrcoef' (Core MATLAB)
% corrcoef returns a matrix, we don't need a toolbox for this
C = corrcoef(Ranks); 

% 5. Manual PRCC Calculation using Inverse Matrix
% Formula: PRCC(i,y) = -C_inv(i,y) / sqrt(C_inv(i,i) * C_inv(y,y))
C_inv = inv(C);
rho = zeros(1, N_vars-1);

for j = 1:N_vars-1
    rho(j) = -C_inv(j, N_vars) / sqrt(C_inv(j,j) * C_inv(N_vars, N_vars));
end

% 6. Create PRCC Results Table
PRCC_Results = table(paramNames', rho', ...
    'VariableNames', {'Parameter', 'PRCC_Coefficient'});

disp('--- PRCC Results (Pure Manual) ---')
disp(PRCC_Results)

%% Updated Plotting Section for Manual Rho
figure('Name', 'GSA Sensitivity Analysis');
subplot(2,1,1);
bar(Sensitivity);
set(gca, 'XTick', 1:length(paramNames), 'XTickLabel', paramNames);
xtickangle(45);
ylabel('Correlation Coefficient');
title('Simple Correlation with Viral Titer');
grid on;

subplot(2,1,2);
bar(rho);
set(gca, 'XTick', 1:length(paramNames), 'XTickLabel', paramNames);
xtickangle(45);
ylabel('PRCC Coefficient');
title('PRCC Analysis (Manual Calculation)');
grid on;
