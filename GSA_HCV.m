clear; close all; clc;

%% 1. Load LHS samples
LHS_table = readtable('LHS_HCV_parameters.csv');
N = height(LHS_table);

%% 2. Load baseline parameter structure
load('param_HCV.mat') 

%% 3. Load steady state initial conditions
load('SteadyState_120h.mat','Yss')
y0 = Yss(end,:);
y0(2) = 1; % Start infection (Initial virus particle)

%% 4. Simulation settings
I_n = 1; I_a = 1; VC = 0.58;
tend = 96*60; % 96 hours in minutes
tspan = linspace(0,tend,1000);

%% 5. Storage
VT_96h = zeros(N,1);
nAOC   = zeros(N,1);
Feasible = zeros(N,1);

%% 6. Thresholds
VT_threshold = 100;
nAOC_threshold = 50;

%% 7. GSA loop
% 25k simulations kabatti 'disp' interval marchithe screen clutter avvadu
for i = 1:N
    if mod(i,100) == 0 || i == 1
        fprintf('Processing Simulation: %d / %d (%.2f%%)\n', i, N, (i/N)*100);
    end

    param_i = param; % Copy baseline structure

    %% Update ALL 20 parameters from LHS table (Correct Mapping)
    % Viral Kinetics
    param_i.k_t            = LHS_table.ktV(i);
    param_i.k_c            = LHS_table.kcV(i);
    param_i.k_r            = LHS_table.krV(i);
    param_i.NCV            = LHS_table.NCV(i);
    param_i.tau            = LHS_table.tau(i);
    param_i.k_en           = LHS_table.kenV(i);
    param_i.k_a            = LHS_table.kaV(i);
    param_i.mu_p           = LHS_table.mu_pV(i);

    % Immune Response / Signaling
    param_i.k_MAVS         = LHS_table.kMAVS(i);
    param_i.kinhISGn       = LHS_table.kinhISGn(i);
    param_i.muISGRNA       = LHS_table.muISGRNA(i);
    param_i.ktISGRNA       = LHS_table.ktISGRNA(i);
    param_i.k_RIGI         = LHS_table.kRIGI(i);
    param_i.k72            = LHS_table.k72(i);
    param_i.degARCBySOCS   = LHS_table.degARCISGn(i);
    param_i.gamma_RIGI     = LHS_table.gammaRIGI(i);

    % Metabolic / Resource Parameters (Varying as per your new LHS)
    param_i.M              = LHS_table.M(i);
    param_i.Omega          = LHS_table.Omega(i);
    param_i.threshold      = LHS_table.threshold(i);
    param_i.n              = LHS_table.nHill(i);

    %% 8. Solve ODE
    try
        % ode23s is good for stiff viral models
        [T,Y] = ode23s(@(t,y) ODEs(t,y,param_i,I_n,I_a,VC), tspan, y0);

        VT_val = Y(end,1); % Viral Titer at 96h
        VT_96h(i) = VT_val;

        % Area Under Curve (nAOC) calculation
        nAOC(i) = trapz(T,Y(:,1)) / (T(end)-T(1));

        % Feasibility check
        if VT_96h(i) >= VT_threshold && nAOC(i) >= nAOC_threshold
            Feasible(i) = 1;
        end
    catch
        % If simulation fails (Singular matrix or timeout)
        VT_96h(i) = NaN;
        nAOC(i) = NaN;
        Feasible(i) = 0;
    end
end

%% 9. Save results table
Results = LHS_table;
Results.VT_96h = VT_96h;
Results.nAOC = nAOC;
Results.Feasible = Feasible;

writetable(Results,'GSA_HCV_Results.csv');
disp('--- Simulation Phase Completed Successfully ---');