clear; close all; clc;

%% Load parameters
loadedData = load('param_JEV.mat');

param = loadedData.param;

%% Add estimated parameters
param.M = 500;
param.Omega = 50;
param.gamma_RIGI = 1;
param.threshold = 500;
param.n = 2;

%% Parameters for LSA
paramNames = {
    'k_t';
    'k_c';
    'k_r';
    'tau';
    'k_en';
    'k_a';
    'k_MAVS';
    'k_RIGI';
    'k72';
    'degARCBySOCS';
    'mu_p';
    'M';
    'Omega';
    'gamma_RIGI'
};

%% Initial conditions
load('SteadyState_120h.mat','Yss');
y0 = Yss(end,:);
y0(2) = 1; % Virus Input

%% Simulation settings
I_n = 1;
I_a = 1;
VC = 1;

tend = 96*60;
tspan = linspace(0,tend,1000);

%% Baseline simulation
[Tbase,Ybase] = ode23s(@(t,y) ODEs(t,y,param,I_n,I_a,VC), tspan, y0);
baseline_output = Ybase(end,1);

%% Storage
Sensitivity = zeros(length(paramNames),1);

%% LSA loop
for i = 1:length(paramNames)
    pname = paramNames{i};
    
    % Ensure the parameter exists in the structure before modifying
    if isfield(param, pname)
        param_new = param;
        oldval = param.(pname);
        
        % +10% perturbation
        param_new.(pname) = 1.1 * oldval;
        
        % Simulate
        [T,Y] = ode23s(@(t,y) ODEs(t,y,param_new,I_n,I_a,VC), tspan, y0);
        new_output = Y(end,1);
        
        % Normalized sensitivity
        if oldval ~= 0
            Sensitivity(i) = ((new_output - baseline_output)/baseline_output) ...
                / ((1.1*oldval - oldval)/oldval);
        else
            Sensitivity(i) = 0; % Avoid division by zero
        end
    else
        warning('Field %s not found in parameter structure.', pname);
    end
end

%% Plot
figure
bar(Sensitivity)
set(gca,'XTick',1:length(paramNames))
set(gca,'XTickLabel',paramNames)
xtickangle(45)

ylabel('Sensitivity Index')
title('Local Sensitivity Analysis')
grid on

%% Display table
Results = table(paramNames, Sensitivity, 'VariableNames', {'Parameter', 'Sensitivity'});
disp(Results)