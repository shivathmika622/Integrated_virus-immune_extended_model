clear; close all; clc;

%% Load parameters
load('param_JEV.mat')

%% Add estimated parameters
param.M = 500;
param.Omega = 50;
param.gamma_RIGI = 1;
param.threshold = 500;
param.n = 2;

%% Initial conditions
load('SteadyState_120h.mat','Yss')

y0 = Yss(end,:);
y0(2) = 1;   % Virus input

%% Simulation settings
I_n = 1;
I_a = 1;
VC  = 1;

tend = 96*60;
tspan = linspace(0,tend,1000);

%% Parameters to sweep
paramList = {'k_r','k_a','M'};

%% Number of sweep points
Npts = 100;

%% Loop over parameters
for p = 1:length(paramList)

    pname = paramList{p};

    %% Baseline value
    oldval = param.(pname);

    %% Sweep range
    sweep_vals = logspace(log10(0.01*oldval), ...
                          log10(10*oldval), ...
                          Npts);

    %% Storage
    VT_96h = zeros(Npts,1);
    nAOC   = zeros(Npts,1);

    %% Sweep loop
    for i = 1:Npts

        param_i = param;

        % Change parameter
        param_i.(pname) = sweep_vals(i);

        %% Solve ODE
        [T,Y] = ode23s(@(t,y) ODEs(t,y,param_i,I_n,I_a,VC), ...
                        tspan, y0);

        %% Viral titer
        VT_96h(i) = Y(end,1);

        %% nAOC
        nAOC(i) = trapz(T,Y(:,1)) / (T(end)-T(1));

    end

    %% Plot VT
    figure

    plot(sweep_vals, VT_96h, '-o', ...
        'LineWidth',1.5, ...
        'MarkerSize',4)

    hold on

    xline(oldval,'--r',['Baseline ', pname]);

    set(gca,'XScale','log','YScale','log')

    xlabel(pname)
    ylabel('V_E at 96 h')

    title(['Parameter Sweep: ', pname, ' vs Viral Titer'])

    grid on

    %% Plot nAOC
    figure

    plot(sweep_vals, nAOC, '-o', ...
        'LineWidth',1.5, ...
        'MarkerSize',4)

    hold on

    xline(oldval,'--r',['Baseline ', pname]);

    set(gca,'XScale','log','YScale','log')

    xlabel(pname)
    ylabel('nAOC')

    title(['Parameter Sweep: ', pname, ' vs nAOC'])

    grid on

end