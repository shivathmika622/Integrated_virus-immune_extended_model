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
warning('off', 'MATLAB:Axes:NegativeDataInLogAxis'); 

% ---------------------------------------------------------
% Parameter Settings (M is Fixed to 1000 here)
% ---------------------------------------------------------
default_Omega = 0;      
default_threshold = 10; 
hill_n = 2; 
M_limit = 1000; % <-- M Value Fixed at 1000

% Parameter Array: [param(1-129), Omega(130), threshold(131), n(132), M(133)]

param.Omega = default_Omega;
param.threshold = default_threshold;
param.n = hill_n;
param.M = M_limit;

param_ss = param;  
tspan_ss = linspace(0, 120 * 60);
[Tss, Yss] = ode23s(@(t, y) ODEs(t, y, param_ss, 1, 0, 0), tspan_ss, Init_Cond);

y0_master = Yss(end,:);
y0_master(2) = 1; % Introduce Virus Input
tspan_inf = linspace(0, 96*60); 

cases = [0 0 0; 1 0 0; 1 1 0; 1 1 1];
case_titles = {'Case 1: No Immune', 'Case 2: Immune Only', 'Case 3: Regulated', 'Case 4: Fully Interacting'};
custom_colors = [0 0 1; 1 0.08 0.58; 0 0.5 0; 1 0.5 0; 0.5 0 0.5; 0 0.75 0.75; 1 0 0];

%% ========================================================================
%% ANALYSIS 1: SWEEPING OMEGA (Silencing Rate)
%% ========================================================================
disp('==================================================================');
disp('   ANALYSIS 1: OMEGA SWEEP (P_NS MAX VALUES)');
disp('==================================================================');
omega_vals = [0, 0.1, 1, 10, 50];
fixed_thresh = 10;

f1 = figure(101); set(f1, 'Name', 'Omega Sweep - ExtVirus', 'units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.9]);
f2 = figure(102); set(f2, 'Name', 'Omega Sweep - ISGav', 'units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.9]);
f3 = figure(103); set(f3, 'Name', 'Omega Sweep - B_U (Open Sites)', 'units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.9]);

for c = 1:size(cases,1)
    I_a = cases(c,1); I_n = cases(c,2); VC  = cases(c,3); 
    fprintf('\n--- %s ---\n', case_titles{c});
    
    % Baseline (Omega = 0, no silencing)
    param.Omega = 0;
    param.threshold = fixed_thresh;
    param.n = hill_n;
    param.M = M_limit;
    [T_base, Y_base] = ode23s(@(t,y) ODEs(t, y, param, I_n, I_a, VC), tspan_inf, y0_master);
    fprintf('Baseline (Omega=0) -> Max P_NS: %.4f nM\n', max(Y_base(:,7)));
    
    figure(101); subplot(2, 2, c); hold on; plot(T_base/60, Y_base(:,1), 'k:', 'LineWidth', 3, 'DisplayName', 'Baseline (\Omega=0)');
    figure(102); subplot(2, 2, c); hold on; plot(T_base/60, Y_base(:,70), 'k:', 'LineWidth', 3, 'DisplayName', 'Baseline (\Omega=0)');
    figure(103); subplot(2, 2, c); hold on; plot(T_base/60, Y_base(:,63), 'k:', 'LineWidth', 3, 'DisplayName', 'Baseline (\Omega=0)');
    
    for m = 1:length(omega_vals)
        % Baseline (Omega = 0)
        param.Omega = omega_vals(m);
        [T, Y] = ode23s(@(t,y) ODEs(t, y, param, I_n, I_a, VC), tspan_inf, y0_master); 
        
        
        fprintf('Omega = %-5g        -> Max P_NS: %.4f nM\n', omega_vals(m), max(Y(:,7)));
        
        leg_str = sprintf('\\Omega = %.2g', omega_vals(m));
        figure(101); subplot(2, 2, c); plot(T/60, Y(:,1), '-', 'LineWidth', 2, 'Color', custom_colors(m,:), 'DisplayName', leg_str);
        figure(102); subplot(2, 2, c); plot(T/60, Y(:,70), '-', 'LineWidth', 2, 'Color', custom_colors(m,:), 'DisplayName', leg_str);
        figure(103); subplot(2, 2, c); plot(T/60, Y(:,63), '-', 'LineWidth', 2, 'Color', custom_colors(m,:), 'DisplayName', leg_str);
    end
    
    for fig = [101, 102, 103]
        figure(fig); subplot(2, 2, c); set(gca, 'YScale', 'log'); 
        xlabel('Time (hours)', 'FontWeight', 'bold'); 
        title(sprintf('%s [%d %d %d]', case_titles{c}, I_a, I_n, VC), ...
            'FontSize', 12, 'FontWeight', 'bold');
        xlim([4 96]); xticks([4 8 16 24 36 48 60 72 84 96]); grid on; legend('Location', 'eastoutside');
    end
    figure(101); subplot(2,2,c); ylabel('ExtVirus (nM)', 'FontWeight', 'bold');
    figure(102); subplot(2,2,c); ylabel('ISGav (nM)', 'FontWeight', 'bold');
    figure(103); subplot(2,2,c); ylabel('B_U (Binding Sites)', 'FontWeight', 'bold');
end
figure(101); sgtitle('Impact of Silencing Power (\Omega) on Extracellular Virus', 'FontSize', 16, 'FontWeight', 'bold');
figure(102); sgtitle('Impact of Silencing Power (\Omega) on Antiviral ISG', 'FontSize', 16, 'FontWeight', 'bold');
figure(103); sgtitle('Depletion of Open DNA Binding Sites (B_U) by \Omega', 'FontSize', 16, 'FontWeight', 'bold');

%% ========================================================================
%% ANALYSIS 2: SWEEPING THRESHOLD (Viral Sensitivity)
%% ========================================================================
disp(' ');
disp('==================================================================');
disp('   ANALYSIS 2: THRESHOLD SWEEP (P_NS MAX VALUES)');
disp('==================================================================');

thresh_vals = [0.5, 10, 100, 300, 500];
fixed_omega = 1;

f5 = figure(201); set(f5, 'Name', 'Threshold Sweep - ExtVirus', 'units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.9]);
f6 = figure(202); set(f6, 'Name', 'Threshold Sweep - ISGav', 'units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.9]);
f7 = figure(203); set(f7, 'Name', 'Threshold Sweep - B_U (Open Sites)', 'units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.9]);

for c = 1:size(cases,1)
    I_a = cases(c,1); I_n = cases(c,2); VC  = cases(c,3); 
    fprintf('\n--- %s ---\n', case_titles{c});
    
    % Baseline (Omega = 0, no silencing)
   
    param.Omega = 0;
    param.threshold = 10;
    param.n = hill_n;
    param.M = M_limit;
    [T_base, Y_base] = ode23s(@(t,y) ODEs(t, y, param, I_n, I_a, VC), tspan_inf, y0_master);
    fprintf('Baseline (No Silence) -> Max P_NS: %.4f nM\n', max(Y_base(:,7)));
    
    figure(201); subplot(2, 2, c); hold on; plot(T_base/60, Y_base(:,1), 'k:', 'LineWidth', 3, 'DisplayName', 'Baseline (\Omega=0)');
    figure(202); subplot(2, 2, c); hold on; plot(T_base/60, Y_base(:,70), 'k:', 'LineWidth', 3, 'DisplayName', 'Baseline (\Omega=0)');
    figure(203); subplot(2, 2, c); hold on; plot(T_base/60, Y_base(:,63), 'k:', 'LineWidth', 3, 'DisplayName', 'Baseline (\Omega=0)');
    
    for m = 1:length(thresh_vals)
       
        param.Omega = fixed_omega;
        param.threshold = thresh_vals(m);
        
        [T, Y] = ode23s(@(t,y) ODEs(t, y, param, I_n, I_a, VC), tspan_inf, y0_master);
        
        fprintf('Thresh = %-5g        -> Max P_NS: %.4f nM\n', thresh_vals(m), max(Y(:,7)));
        
        leg_str = sprintf('Thresh = %g', thresh_vals(m));
        figure(201); subplot(2, 2, c); hold on; plot(T/60, Y(:,1), '-', 'LineWidth', 2, 'Color', custom_colors(m,:), 'DisplayName', leg_str);
        figure(202); subplot(2, 2, c); hold on; plot(T/60, Y(:,70), '-', 'LineWidth', 2, 'Color', custom_colors(m,:), 'DisplayName', leg_str);
        figure(203); subplot(2, 2, c); hold on; plot(T/60, Y(:,63), '-', 'LineWidth', 2, 'Color', custom_colors(m,:), 'DisplayName', leg_str);
    end
    
    for fig = [201, 202, 203]
        figure(fig); subplot(2, 2, c); set(gca, 'YScale', 'log'); 
        xlabel('Time (hours)', 'FontWeight', 'bold'); 
        title(sprintf('%s [%d %d %d]', case_titles{c}, I_a, I_n, VC), ...
            'FontSize', 12, 'FontWeight', 'bold');
        xlim([4 96]); xticks([4 8 16 24 36 48 60 72 84 96]); grid on; legend('Location', 'eastoutside');
    end
    figure(201); subplot(2,2,c); ylabel('ExtVirus (nM)', 'FontWeight', 'bold');
    figure(202); subplot(2,2,c); ylabel('ISGav (nM)', 'FontWeight', 'bold');
    figure(203); subplot(2,2,c); ylabel('B_U (Binding Sites)', 'FontWeight', 'bold');
end
figure(201); sgtitle('Impact of Silencing Threshold on Virus', 'FontSize', 16, 'FontWeight', 'bold');
figure(202); sgtitle('Impact of Silencing Threshold on Antiviral ISG', 'FontSize', 16, 'FontWeight', 'bold');
figure(203); sgtitle('Impact of Silencing Threshold on Binding Sites (B_U)', 'FontSize', 16, 'FontWeight', 'bold');

disp('✅ Simulation Complete! All 6 Figures generated successfully.');
