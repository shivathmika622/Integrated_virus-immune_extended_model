clear; close all; clc;

%% Load GSA results
Results = readtable('GSA_JEV_Results.csv');

%% Parameters to analyze
paramNames = {
    'ktV'
    'kcV'
    'krV'
    'tau'
    'kenV'
    'kaV'
    'kMAVS'
    'kRIGI'
    'k72'
    'degARCISGn'
    'mu_pV'
    'M'
    'Omega'
    'gammaRIGI'
    };

%% Storage
Sensitivity = zeros(length(paramNames),1);

%% Correlation-based GSA
for i = 1:length(paramNames)

    x = Results.(paramNames{i});
    y = Results.VT_96h;

    R = corrcoef(x,y,'Rows','complete');

    Sensitivity(i) = R(1,2);

end

%% Plot
figure

bar(Sensitivity)

set(gca,'XTick',1:length(paramNames))
set(gca,'XTickLabel',paramNames)

xtickangle(45)

ylabel('Correlation with VT\_96h')

title('Global Sensitivity Analysis')

grid on

%% Display table
GSATable = table(paramNames, Sensitivity, ...
    'VariableNames', {'Parameter','Sensitivity'});

disp(GSATable)

%% Save
writetable(GSATable,'GSA_Sensitivity_Table.csv')

disp('GSA analysis completed')