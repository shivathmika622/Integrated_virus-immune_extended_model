clear; close all; clc;

%% ========================================================================
%% SECTION 1: INITIAL CONDITIONS & SETUP
%% ========================================================================

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

load('param_HCV.mat')

% Steady state run with default parameters (gamma_RIGI = 0)
param_base = param;
param_base.gamma_RIGI = 0;  
tspan_ss = linspace(0, 120 * 60);

[Tss, Yss] = ode23s(@(t, y) ODEs(t, y, param_base, 1, 0, 0), tspan_ss, Init_Cond);

y0_master = Yss(end,:);
y0_master(2) = 1; % Introduce Virus Input
tspan_inf = linspace(0, 96*60);

%% ========================================================================
%% SECTION 2: RIG-I SUPPRESSION SWEEP (gamma_RIGI = 0.005 to 1)
%% ========================================================================

% Sweep range — 5 values only (removed 5 and 10)
gamma_vals = [0.005, 0.05, 0.1, 0.5, 1];

% Corresponding colors — 5 colors only (removed Cyan and Red)
custom_colors = [
    0      0      1   ;   % 0.005 → Blue
    1      0.08   0.58;   % 0.05  → Hot Pink
    0      0.5    0   ;   % 0.1   → Green
    1      0.5    0   ;   % 0.5   → Orange
    0.5    0      0.5 ;   % 1     → Purple
    ];

cases = [0 0 0; 1 0 0; 1 1 0; 1 1 1];
case_titles = {'Case 1: No Immune [0 0 0]', 'Case 2: Immune Only [1 0 0]', ...
    'Case 3: Regulated [1 1 0]', 'Case 4: Fully Interacting [1 1 1]'};

f_sw1 = figure(101); set(f_sw1, 'Name', 'ExtVirus Analysis', 'units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.9]);
f_sw2 = figure(102); set(f_sw2, 'Name', 'ISGav Analysis', 'units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.9]);

disp('===================================================');
disp('   P_NS (NON-STRUCTURAL PROTEIN) MAX LEVELS        ');
disp('===================================================');

for c = 1:size(cases,1)
    I_a = cases(c,1); I_n = cases(c,2); VC  = cases(c,3);
    fprintf('\n--- %s ---\n', case_titles{c});

    % Baseline Dotted Line (No Suppression, gamma = 0)
    param_unlim = param;
    param_unlim.gamma_RIGI = 0;
    [T_base, Y_base] = ode23s(@(t,y) ODEs(t, y, param_unlim, I_n, I_a, VC), tspan_inf, y0_master);

    fprintf('Baseline (No Suppression) -> Max P_NS: %.4f nM\n', max(Y_base(:,7)));

    figure(101); subplot(2, 2, c); hold on; plot(T_base/60, Y_base(:,1), 'k:', 'LineWidth', 3, 'DisplayName', 'Baseline (\gamma=0)');
    figure(102); subplot(2, 2, c); hold on; plot(T_base/60, Y_base(:,70), 'k:', 'LineWidth', 3, 'DisplayName', 'Baseline (\gamma=0)');

    % Loop through the gamma values
    for m = 1:length(gamma_vals)
            param_sweep = param;
            param_sweep.gamma_RIGI = gamma_vals(m);

            [T, Y] = ode23s(@(t,y) ODEs(t, y, param_sweep, I_n, I_a, VC), tspan_inf, y0_master);

        fprintf('Gamma_RIGI = %-5g       -> Max P_NS: %.4f nM\n', gamma_vals(m), max(Y(:,7)));

        legend_str = sprintf('\\gamma_{RIGI} = %.3g', gamma_vals(m));

        figure(101); subplot(2, 2, c); plot(T/60, Y(:,1), '-', 'LineWidth', 2, 'Color', custom_colors(m,:), 'DisplayName', legend_str);
        figure(102); subplot(2, 2, c); plot(T/60, Y(:,70), '-', 'LineWidth', 2, 'Color', custom_colors(m,:), 'DisplayName', legend_str);
    end

    % Format Subplots
    for fig_num = 101:102
        figure(fig_num); subplot(2, 2, c);
        set(gca, 'YScale', 'log'); xlabel('Time (hours)', 'FontWeight', 'bold');
        title(case_titles{c}, 'FontSize', 12); xlim([4 96]); xticks([4 8 16 24 36 48 60 72 84 96]);
        grid on; legend('Location', 'eastoutside');
    end
    figure(101); subplot(2,2,c); ylabel('ExtVirus (nM)', 'FontWeight', 'bold');
    figure(102); subplot(2,2,c); ylabel('ISGav (nM)', 'FontWeight', 'bold');
end

figure(101); sgtitle('Impact of RIG-I Suppression (\gamma_{RIGI}) on Extracellular Virus', 'FontSize', 16, 'FontWeight', 'bold');
figure(102); sgtitle('Impact of RIG-I Suppression (\gamma_{RIGI}) on Antiviral ISG', 'FontSize', 16, 'FontWeight', 'bold');
