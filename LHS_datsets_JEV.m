clc;
clear;

%% Number of samples
N = 25000;

%% Nominal parameter values (JEV + estimated values)
params.ktV           = 1.6e2;        % kt,V
params.kcV           = 1.6e-2;       % kc,V
params.krV           = 3.7;          % kr,V
params.NCV           = 1.21e3;       % NC,V
params.tau           = 2.6;          % tau
params.kenV          = 0.43;         % ken,V
params.kaV           = 8e-9;         % ka,V
params.kMAVS         = 9e-3;         % kMAVS
params.kinhISGn      = 889.4;        % kinhISGn
params.muISGRNA      = 2.8e-3;       % muISGRNA
params.ktISGRNA      = 1.2e-4;       % kt,ISGRNA
params.M             = 500;          % estimated
params.Omega         = 50;           % estimated
params.gammaRIGI     = 1;            % estimated
params.kRIGI         = 0.01;         % kRIGI
params.k72           = 0.147;        % k72
params.degARCISGn    = 0.0147;       % degARCISGn
params.mu_pV         = 0.11;         % mup,V

params.threshold     = 500;          % estimated
params.nHill         = 2;            % estimated

%% Parameter names
paramNames = fieldnames(params);
paramNames = paramNames(:)'; % Ensure it is a row cell array for table creation

%% Convert to vector (enforce row vector orientation)
p0 = struct2array(params);
p0 = p0(:)';

%% Number of parameters
numParams = length(p0);

%% LHS sampling
lhs = zeros(N, numParams);

for j = 1:numParams

    % Divide interval into N bins
    temp = ((0:N-1)' + rand(N,1)) / N;

    % Random permutation
    lhs(:,j) = temp(randperm(N));

end

%% Lower and upper bounds (1 order of magnitude)
lb = log10(p0 / 10);
ub = log10(p0 * 10);

%% Scale samples (log-uniform)
samples_log = lb + lhs .* (ub - lb);

%% Convert back to linear space
samples = 10.^samples_log;

%% Convert to table
LHS_table = array2table(samples, ...
    'VariableNames', paramNames);

%% Save
writetable(LHS_table, 'LHS_JEV_parameters.csv');

disp('LHS sampling completed and saved to LHS_JEV_parameters.csv');