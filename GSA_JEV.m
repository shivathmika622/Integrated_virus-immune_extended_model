clear; close all; clc;

%% Load LHS samples
LHS_table = readtable('LHS_JEV_parameters.csv');

%% Number of parameter sets
N = height(LHS_table);

%% Load baseline parameter structure
load('param_JEV.mat')

%% Add fixed estimated parameters (If not already in table)
param.M = 500;
param.Omega = 50;
param.gamma_RIGI = 1;
param.threshold = 500;
param.n = 2;

%% Load steady state initial conditions
load('SteadyState_120h.mat','Yss')
y0 = Yss(end,:);
y0(2) = 1; % Start infection

%% Simulation settings
I_n = 1; I_a = 1; VC = 1;
tend = 96*60;
tspan = linspace(0,tend,1000);

%% Storage
VT_96h = zeros(N,1);
nAOC   = zeros(N,1);
Feasible = zeros(N,1);

%% Thresholds
VT_threshold = 100;
nAOC_threshold = 50;

%% GSA loop
for i = 1:N
    disp(['Simulation: ', num2str(i), '/', num2str(N)])

    param_i = param;

    %% Update ALL 20 parameters from LHS table
    param_i.k_t            = LHS_table.ktV(i);
    param_i.k_c            = LHS_table.kcV(i);
    param_i.k_r            = LHS_table.krV(i);
    param_i.NCV            = LHS_table.NCV(i);       % Added
    param_i.tau            = LHS_table.tau(i);
    param_i.k_en           = LHS_table.kenV(i);
    param_i.k_a            = LHS_table.kaV(i);
    param_i.k_MAVS         = LHS_table.kMAVS(i);
    param_i.kinhISGn       = LHS_table.kinhISGn(i);  % Added
    param_i.muISGRNA       = LHS_table.muISGRNA(i);  % Added
    param_i.ktISGRNA       = LHS_table.ktISGRNA(i);  % Added
    param_i.k_RIGI         = LHS_table.kRIGI(i);
    param_i.k72            = LHS_table.k72(i);
    param_i.degARCBySOCS   = LHS_table.degARCISGn(i);
    param_i.mu_p           = LHS_table.mu_pV(i);
    param_i.M              = LHS_table.M(i);
    param_i.Omega          = LHS_table.Omega(i);
    param_i.gamma_RIGI     = LHS_table.gammaRIGI(i);
    param_i.threshold      = LHS_table.threshold(i);
    param_i.n              = LHS_table.nHill(i);

    %% Solve ODE
    try
        [T,Y] = ode23s(@(t,y) ODEs(t,y,param_i,I_n,I_a,VC), tspan, y0);
        VT_96h(i) = Y(end,1);
        nAOC(i) = trapz(T,Y(:,1)) / (T(end)-T(1));

        if VT_96h(i) >= VT_threshold && nAOC(i) >= nAOC_threshold
            Feasible(i) = 1;
        end
    catch
        VT_96h(i) = NaN;
        nAOC(i) = NaN;
        Feasible(i) = 0;
    end
end

%% Save results table
Results = LHS_table;
Results.VT_96h = VT_96h;
Results.nAOC = nAOC;
Results.Feasible = Feasible;
writetable(Results,'GSA_JEV_Results.csv')
disp('Simulation Phase Completed.')
