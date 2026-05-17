clear; close all; clc;

%% 1. Initial Conditions Vector
V_E = 0; V_0 = 0; V_I = 0; R_cyt = 0; R_CM = 0; P_S = 0; P_NS = 0; RC_CM = 0; R_ds = 0;
RIGI = 5.3422; aRIGI = 0; MAVS = 277.78; aMAVS = 0; IKKe = 3.0832; pIKKe = 0; TBK1 = 97.169; pTBK1 = 0;
IRF3 = 37.862; pIRF3 = 0; IKK = 37.969; aIKK = 0; NFkB_IkBac = 11.3587; pNFkBn = 0; NFkBn = 0;
NFkBc = 101.735; IkBac = 0; IRF7 = 24.92; pIRF7 = 0; IFNb_mRNA = 0; IFNa_mRNA = 0; IFNl_mRNA = 0;
IFN_c = 0; IFNl_c = 0; JAK = 151.88; RJC = 0; STAT1c = 1114.8; CP = 20; ISGn = 0; IFNex = 0; 
STAT2c = 6.5019; TYK = 20.701; RTC = 0; ARC = 0; IFNAR1 = 1000; IFNAR2 = 1000; IFNAR_d = 0;
IRF9c = 45; ARC_STAT2c = 0; ARC_STAT12c = 0; STAT2c_IRF9 = 0; ISGF3c = 0; PSC_c = 0; ISGF3_CP = 0;
PSC_CP = 0; NP = 40; STAT1n = 0; STAT2n = 0; PIAS = 41.96; PSC_n = 0; IRF9n = 0; ISGF3n = 0;
PSC_NP = 0; B_u = 500; B_o_NP = 0; B_o = 0; ISGF3_PIAS = 0; STAT2n_IRF9 = 0; ISGF3_NP = 0;
ISGav_mRNA = 0; ISGav = 0; ISGn_mRNA_n = 0; IRF9_mRNA_n = 0; IRF7_mRNA = 0; ISGn_mRNA_c = 0; IRF9_mRNA_c = 0;

Init_Cond = [V_E V_0 V_I R_cyt R_CM P_S P_NS RC_CM R_ds RIGI aRIGI MAVS aMAVS IKKe pIKKe TBK1 ...
    pTBK1 IRF3 pIRF3 IKK aIKK NFkB_IkBac pNFkBn NFkBn NFkBc IkBac IRF7 pIRF7 IFNb_mRNA IFNa_mRNA IFNl_mRNA IFN_c IFNl_c JAK RJC STAT1c CP ISGn IFNex STAT2c TYK RTC ARC IFNAR1 IFNAR2 IFNAR_d IRF9c ARC_STAT2c ARC_STAT12c STAT2c_IRF9 ISGF3c PSC_c ISGF3_CP PSC_CP NP STAT1n STAT2n PIAS PSC_n IRF9n ISGF3n PSC_NP B_u B_o_NP B_o ISGF3_PIAS STAT2n_IRF9 ISGF3_NP ISGav_mRNA ISGav ISGn_mRNA_n IRF9_mRNA_n IRF7_mRNA ISGn_mRNA_c IRF9_mRNA_c];

%% 2. Setup Parameters & Variables
load('param_HCV.mat') 

M_vals = [0.005, 0.5, 1, 10, 100, 1000, 5000];

custom_colors = [
    0, 0, 1;          % Blue
    1, 0.08, 0.58;    % Hot Pink
    0, 0.5, 0;        % Dark Green
    1, 0.5, 0;        % Orange
    0.5, 0, 0.5;      % Purple
    0, 0.75, 0.75;    % Cyan
    1, 0, 0           % Red
    ];

tspan_ss = linspace(0, 120*60);
tspan = linspace(0, 96*60); 

cases = [0 0 0; 1 0 0; 1 1 0; 1 1 1];
case_titles = {'Case 1: No Immune [0 0 0]', 'Case 2: Immune Only [1 0 0]', ...
    'Case 3: Regulated [1 1 0]', 'Case 4: Fully Interacting [1 1 1]'};

%% 3. Establish Base Steady State (RUNS EXACTLY ONCE)

param.M = 1000;            
param_ss = param;           
[~, Yss] = ode23s(@(t,y) ODEs(t, y, param_ss, 1, 0, 0), tspan_ss, Init_Cond);

y0_master = Yss(end, :);
y0_master(2) = 1; % Introduce Virus Input

%% 4. Setup Figures (Exactly 8 Graphs across 2 Windows)
f1 = figure('Name', 'ExtVirus Analysis', 'units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.9]);
f2 = figure('Name', 'ISGav Analysis', 'units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.9]);

disp('---------------------------------------------------');
disp('Calculated P_NS (Protein) Max Values:');
disp('---------------------------------------------------');

%% 5. Master Loop: Simulate and Plot everything
for c = 1:size(cases,1)

    I_a = cases(c,1); I_n = cases(c,2); VC  = cases(c,3);

    fprintf('\n--- %s ---\n', case_titles{c});

    % Baseline Dotted Line
    
    param.M = 1e6;                  
    baseline_param = param;
    [T_base, Y_base] = ode23s(@(t,y) ODEs(t, y, baseline_param, I_n, I_a, VC), tspan, y0_master);

    fprintf('Baseline (No Limit) -> Max P_NS: %.4f nM\n', max(Y_base(:,7)));

    figure(f1); subplot(2, 2, c); hold on; plot(T_base/60, Y_base(:,1), 'k:', 'LineWidth', 3, 'DisplayName', 'Baseline (No M Limit)');
    figure(f2); subplot(2, 2, c); hold on; plot(T_base/60, Y_base(:,70), 'k:', 'LineWidth', 3, 'DisplayName', 'Baseline (No M Limit)');

    % Loop through limited M values
    for m = 1:length(M_vals)
      
        param.M = M_vals(m);                  
        current_param = param; 
        [T, Y] = ode23s(@(t,y) ODEs(t, y, current_param, I_n, I_a, VC), tspan, y0_master);

        fprintf('M = %g -> Max P_NS: %.4f nM\n', M_vals(m), max(Y(:,7)));

        legend_str = sprintf('M = %.3g', M_vals(m));

        figure(f1); subplot(2, 2, c); plot(T/60, Y(:,1), '-', 'LineWidth', 2, 'Color', custom_colors(m,:), 'DisplayName', legend_str);
        figure(f2); subplot(2, 2, c); plot(T/60, Y(:,70), '-', 'LineWidth', 2, 'Color', custom_colors(m,:), 'DisplayName', legend_str);
    end

    % Format Subplots
    figure(f1); subplot(2, 2, c);
    set(gca, 'YScale', 'log'); ylabel('ExtVirus (nM)', 'FontWeight', 'bold'); xlabel('Time (hours)', 'FontWeight', 'bold');
    title(case_titles{c}, 'FontSize', 12); xlim([4 96]); xticks([4 8 16 24 36 48 60 72 84 96]); grid on; 
    legend('Location', 'eastoutside'); % <--- CHANGED THIS TO 'eastoutside'

    figure(f2); subplot(2, 2, c);
    set(gca, 'YScale', 'log'); ylabel('ISGav (nM)', 'FontWeight', 'bold'); xlabel('Time (hours)', 'FontWeight', 'bold');
    title(case_titles{c}, 'FontSize', 12); xlim([4 96]); xticks([4 8 16 24 36 48 60 72 84 96]); grid on; 
    legend('Location', 'eastoutside'); % <--- CHANGED THIS TO 'eastoutside'
    disp('---------------------------------------------------');
end
