clc;
clear;

%% Number of samples
N = 1000;

%% Nominal parameter values for PV (Poliovirus - From Table 1 & Shared)
params.ktV           = 18.9;         % kt,V (PV Table 1)
params.kcV           = 1.2e-2;       % kc,V (PV Table 1)
params.krV           = 220;          % kr,V (PV Table 1 - Fast replication)
params.NCV           = 1.62e3;       % NC,V (PV Table 1)
params.tau           = 4.2;          % tau  (PV Table 1)
params.kenV          = 0.43;          % ke,V  (PV Table 1 - Export rate)
params.kaV           = 5.6e-3;       % ka,V  (PV Table 1)
params.kMAVS         = 9e-3;         % Shared immune param
params.kinhISGn      = 889.4;        % Shared immune param
params.muISGRNA      = 2.8e-3;       % Shared immune param
params.ktISGRNA      = 1.2e-4;       % Shared immune param
params.M             = 500;          % estimated
params.Omega         = 50;           % estimated
params.gammaRIGI     = 1;            % estimated
params.kRIGI         = 0.01;         % Shared immune param
params.k72           = 0.147;        % Shared immune param
params.degARCISGn    = 0.0147;       % Shared immune param
params.mu_pV         = 0.11;         % Shared immune param
params.threshold     = 500;          % estimated
params.nHill         = 2;            % estimated

%% Parameter names
paramNames = fieldnames(params);
paramNames = paramNames(:)'; 

%% Convert to vector
p0 = struct2array(params);
p0 = p0(:)';

%% Number of parameters
numParams = length(p0);

%% LHS sampling
lhs = zeros(N, numParams);
for j = 1:numParams
    temp = ((0:N-1)' + rand(N,1)) / N;
    lhs(:,j) = temp(randperm(N));
end

%% Lower and upper bounds (1 order of magnitude)
lb = log10(p0 * 0.1);
ub = log10(p0 * 10);

%% Scale samples (log-uniform)
samples_log = lb + lhs .* (ub - lb);

%% Convert back to linear space
samples = 10.^samples_log;

%% Convert to table
LHS_table = array2table(samples, 'VariableNames', paramNames);

%% Save
writetable(LHS_table, 'LHS_PV_parameters.csv');
disp('LHS sampling for Polio Virus (PV) completed and saved to LHS_PV_parameters.csv');