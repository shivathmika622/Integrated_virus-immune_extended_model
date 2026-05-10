clear; close all; clc;

%% 1. Load baseline parameters
loadedData = load('param_HCV.mat');
param = loadedData.param;

%% 2. Add fixed/baseline values for all 20 parameters
param.M           = 500;
param.Omega       = 50;
param.gamma_RIGI  = 1;
param.threshold   = 500;
param.n           = 2;

param.NCV         = 88;
param.kinhISGn    = 889.4;
param.muISGRNA    = 2.8e-3;
param.ktISGRNA    = 1.2e-4;

%% 3. Parameters list for LSA
paramNames = {
    'k_t'; 'k_c'; 'k_r'; 'NCV'; 'tau'; ...
    'k_en'; 'k_a'; 'mu_p'; 'k_MAVS'; 'kinhISGn'; ...
    'muISGRNA'; 'ktISGRNA'; 'k_RIGI'; 'k72'; 'degARCBySOCS'; ...
    'gamma_RIGI'; 'M'; 'Omega'; 'threshold'; 'n'
};

%% 4. Initial conditions
load('SteadyState_120h.mat','Yss');
y0 = Yss(end,:);
y0(2) = 1;   % Start infection

%% 5. Simulation settings
I_n = 1;
I_a = 1;
VC   = 0.58;

tend  = 96 * 60;
tspan = linspace(0, tend, 1000);

%% 6. Baseline simulation
fprintf('Running baseline simulation...\n');
try
    [~, Ybase] = ode23s(@(t,y) ODEs(t,y,param,I_n,I_a,VC), tspan, y0);
    baseline_output = Ybase(end,1);
catch ME
    error('Baseline simulation failed: %s', ME.message);
end

% Avoid division by very small numbers
if abs(baseline_output) < 1e-8
    warning('Baseline output is very small. Sensitivity values may be unstable.');
    baseline_output = 1e-8;
end

%% 7. Central-difference LSA
N_p = length(paramNames);
Sensitivity = nan(N_p,1);

delta = 0.10;   % 10% perturbation

for i = 1:N_p
    pname = paramNames{i};
    fprintf('Analyzing Parameter (%d/%d): %s\n', i, N_p, pname);

    if ~isfield(param, pname)
        warning('Field %s not found in param structure!', pname);
        continue;
    end

    p0 = param.(pname);

    % Skip parameters that cannot be perturbed multiplicatively
    if p0 == 0
        warning('Parameter %s is zero, cannot use multiplicative perturbation.', pname);
        continue;
    end

    % Create plus and minus perturbed parameter sets
    param_plus  = param;
    param_minus = param;

    param_plus.(pname)  = p0 * (1 + delta);
    param_minus.(pname) = p0 * (1 - delta);

    try
        [~, Yp] = ode23s(@(t,y) ODEs(t,y,param_plus,I_n,I_a,VC), tspan, y0);
        [~, Ym] = ode23s(@(t,y) ODEs(t,y,param_minus,I_n,I_a,VC), tspan, y0);

        Y_plus  = Yp(end,1);
        Y_minus = Ym(end,1);

        % Central-difference normalized sensitivity:
        % S = ((Y+ - Y-) / (2*Y0)) / delta
        Sensitivity(i) = ((Y_plus - Y_minus) / (2 * baseline_output)) / delta;

    catch ME
        warning('Simulation failed for %s: %s', pname, ME.message);
        Sensitivity(i) = NaN;
    end
end

%% 8. Plot results
%% 8. Plot results
figure('Position', [100, 100, 1200, 550]);   % wider figure

b = bar(Sensitivity);

ax = gca;
ax.XTick = 1:N_p;
ax.XTickLabel = paramNames;
ax.TickLabelInterpreter = 'none';   % important for underscores / full text
ax.FontSize = 9;

xtickangle(45);
ylabel('Normalized Sensitivity Index');
title('Local Sensitivity Analysis (LSA)');
grid on;

% give more room at the bottom for long labels
ax.Position(2) = 0.22;
ax.Position(4) = 0.68;

% Add value labels on bars
%xtips = b.XEndPoints;
%ytips = b.YEndPoints;
%labels = string(round(Sensitivity, 4));
%text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
   % 'VerticalAlignment', 'bottom', 'FontSize', 8);

%% 9. Display results table
ResultsTable = table(paramNames, Sensitivity, ...
    'VariableNames', {'Parameter', 'SensitivityIndex'});
disp(ResultsTable);