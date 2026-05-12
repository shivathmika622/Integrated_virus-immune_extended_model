% 1. Base template ga HCV file ni load cheyandi
load('param_HCV.mat'); 
param_PV = param; % 134 parameters unna structure copy avthundi

% 2. Ippudu mee daggara unna 20 PV parameters ni override cheyandi
% (Mapping names to your structure)
param_PV.k_t = 18.9;           % ktV
param_PV.k_c = 1.2e-2;         % kcV
param_PV.k_r = 220;            % krV
param_PV.rcsat = 1.62e3;       % NCV
param_PV.tau = 4.2;            % tau
param_PV.k_e = 0.43;           % kenV
param_PV.k_a = 5.6e-3;         % kaV
param_PV.k_MAVS = 9e-3;        
param_PV.kinhBySOCS = 889.4;   % kinhISGn
param_PV.mu_ISG_RNA = 2.8e-3;  % muISGRNA
param_PV.k_transISG = 1.2e-4;  % ktISGRNA
param_PV.M = 500;
param_PV.Omega = 50;
param_PV.gamma_RIGI = 1;
param_PV.k_RIGI = 0.01;
param_PV.k72 = 0.147;
param_PV.degARCBySOCS = 0.0147; % degARCISGn
param_PV.mu_p = 0.11;          % mu_pV
param_PV.threshold = 500;
param_PV.n = 2;                % nHill

% 3. Ee updated structure ni PV file ga save cheyandi
param = param_PV;
save('param_PV.mat', 'param');
disp('param_PV.mat file successfully created with 134 parameters!');