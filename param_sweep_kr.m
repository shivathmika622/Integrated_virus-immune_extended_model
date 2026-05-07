clear; close all; clc;

%% Load baseline parameters
load('param_JEV.mat')

% Fixed estimated parameters
param.M = 500;
param.Omega = 50;
param.gamma_RIGI = 1;
param.threshold = 500;
param.n = 2;

%% Initial conditions
load('SteadyState_120h.mat','Yss')
y0 = Yss(end,:);
y0(2) = 1;   % virus input

%% Switches
I_n = 1;
I_a = 1;
VC  = 1;

%% Time settings
tend = 96*60;
tspan = linspace(0,tend,1000);

%% Sweep one parameter: k_r
oldval = param.k_r;                 % baseline value
Npts   = 100;                       % number of sweep points
k_r_vals = logspace(log10(0.01*oldval), log10(10*oldval), Npts);

%% Storage
VT_96h = zeros(Npts,1);
nAOC   = zeros(Npts,1);

%% Sweep loop
for i = 1:Npts
    param_i = param;
    param_i.k_r = k_r_vals(i);

    [T,Y] = ode23s(@(t,y) ODEs(t,y,param_i,I_n,I_a,VC), tspan, y0);

    VT_96h(i) = Y(end,1);
    nAOC(i)   = trapz(T, Y(:,1)) / (T(end)-T(1));
end

%% Plot 1: final viral titer vs k_r
figure
plot(k_r_vals, VT_96h, '-o', 'LineWidth', 1.5, 'MarkerSize', 4)
hold on
xline(oldval, '--r');
set(gca,'XScale','log','YScale','log')
xlabel('k_r')
ylabel('V_E at 96 h')
title('Parameter Sweep: k_r vs Viral Titer')
grid on

%% Plot 2: nAOC vs k_r
figure
plot(k_r_vals, nAOC, '-o', 'LineWidth', 1.5, 'MarkerSize', 4)
hold on
xline(oldval, '--r');
set(gca,'XScale','log','YScale','log')
xlabel('k_r')
ylabel('nAOC')
title('Parameter Sweep: k_r vs nAOC')
grid on

