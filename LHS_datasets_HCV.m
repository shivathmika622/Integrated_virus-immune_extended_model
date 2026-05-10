clc;
clear;

%% Number of samples
N = 1000;

%% Nominal parameter values (HCV + estimated values)
params.ktV           = 23.7;         % kt,V (HCV specific)
params.kcV           = 2.6e-3;       % kc,V (HCV specific)
params.krV           = 3.6;          % kr,V (HCV specific)
params.NCV           = 88;           % NC,V (HCV specific)
params.tau           = 5.8;          % tau (HCV specific)
params.kenV          = 0.43;         % ken,V (Shared value)
params.kaV           = 3.6e-8;       % ka,V (HCV specific)
params.kMAVS         = 9e-3;         % kMAVS (Shared)
params.kinhISGn      = 889.4;        % kinhISGn (Shared/HCV)
params.muISGRNA      = 2.8e-3;       % muISGRNA (Shared)
params.ktISGRNA      = 1.2e-4;       % kt,ISGRNA (Shared)
params.M             = 500;          % estimated
params.Omega         = 50;           % estimated
params.gammaRIGI     = 1;            % estimated
params.kRIGI         = 0.01;         % kRIGI (Shared)
params.k72           = 0.147;        % k72 (Shared/HCV)
params.degARCISGn    = 0.0147;       % degARCISGn (Shared)
params.mu_pV         = 0.11;         % mup,V (Shared/HCV)
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
lb = log10(p0 * 0.1);
ub = log10(p0 * 10);

%% Scale samples (log-uniform)
samples_log = lb + lhs .* (ub - lb);

%% Convert back to linear space
samples = 10.^samples_log;

%% Convert to table
LHS_table = array2table(samples, ...
    'VariableNames', paramNames);

%% Save
writetable(LHS_table, 'LHS_HCV_parameters.csv');

disp('LHS sampling completed and saved to LHS_HCV_parameters.csv');