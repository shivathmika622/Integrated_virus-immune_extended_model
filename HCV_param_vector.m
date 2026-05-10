%% Import data from spreadsheet
clear; clc; close all; 
%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 129);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:DY2";

% Specify column names and types
opts.VariableNames = ["Vc2n", "Vn2c", "a_RC", "b_IRF3", "b_MAVS", "b_RIGI", "b_prot", "deg", "degARCBySOCS", "degRecBySOCS", "k1", "k10", "k11", "k12", "k13", "k14", "k15", "k16", "k17", "k18", "k19", "k2", "k20", "k21", "k22", "k23", "k27", "k29", "k3", "k31", "k32", "k34", "k35", "k36", "k37", "k38", "k39", "k4", "k40", "k41", "k42", "k43", "k44", "k45", "k46", "k47", "k48", "k49", "k5", "k50", "k51", "k52", "k53", "k54", "k56", "k57", "k58", "k59", "k6", "k60", "k61", "k62", "k63", "k64", "k65", "k66", "k67", "k69", "k7", "k70", "k71", "k72", "k73", "k74", "k75", "k76", "k77", "k78", "k79", "k8", "k9", "k_IFN", "k_IKK", "k_IKKe_TBK1", "k_IRF3_IKKe_TBK1", "k_MAVS", "k_RIGI", "k_TFBS_IFNa", "k_TFBS_IFNb", "k_TFBS_IFNl", "k_a", "k_act", "k_c", "k_deph", "k_e", "k_en", "k_expr_IkBa", "k_f", "k_inh_p65", "k_mRNA_IFNb", "k_mRNA_IFNl", "k_mRNA_MX1", "k_mx1_mRNA", "k_r", "k_rigi_synt", "k_t", "k_transISG", "k_transISGn", "k_trans_IFNl", "k_transp_NFkB", "kinhBySOCS", "ks", "mu_IFN", "mu_IFNl", "mu_IkBa", "mu_V_E", "mu_V_I", "mu_mRNA_IFNa", "mu_mRNA_IFNb", "mu_mRNA_IFNl", "mu_mRNA_mx1", "mu_mx1", "mu_p", "mu_r", "mu_rigi", "nSP", "rcsat", "tau", "tau_5"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
tbl = readtable("Param.xlsx", opts, "UseExcel", false);

%% Convert to output type
Vc2n = tbl.Vc2n;
Vn2c = tbl.Vn2c;
a_RC = tbl.a_RC;
b_IRF3 = tbl.b_IRF3;
b_MAVS = tbl.b_MAVS;
b_RIGI = tbl.b_RIGI;
b_prot = tbl.b_prot;
deg = tbl.deg;
degARCBySOCS = tbl.degARCBySOCS;
degRecBySOCS = tbl.degRecBySOCS;
k1 = tbl.k1;
k10 = tbl.k10;
k11 = tbl.k11;
k12 = tbl.k12;
k13 = tbl.k13;
k14 = tbl.k14;
k15 = tbl.k15;
k16 = tbl.k16;
k17 = tbl.k17;
k18 = tbl.k18;
k19 = tbl.k19;
k2 = tbl.k2;
k20 = tbl.k20;
k21 = tbl.k21;
k22 = tbl.k22;
k23 = tbl.k23;
k27 = tbl.k27;
k29 = tbl.k29;
k3 = tbl.k3;
k31 = tbl.k31;
k32 = tbl.k32;
k34 = tbl.k34;
k35 = tbl.k35;
k36 = tbl.k36;
k37 = tbl.k37;
k38 = tbl.k38;
k39 = tbl.k39;
k4 = tbl.k4;
k40 = tbl.k40;
k41 = tbl.k41;
k42 = tbl.k42;
k43 = tbl.k43;
k44 = tbl.k44;
k45 = tbl.k45;
k46 = tbl.k46;
k47 = tbl.k47;
k48 = tbl.k48;
k49 = tbl.k49;
k5 = tbl.k5;
k50 = tbl.k50;
k51 = tbl.k51;
k52 = tbl.k52;
k53 = tbl.k53;
k54 = tbl.k54;
k56 = tbl.k56;
k57 = tbl.k57;
k58 = tbl.k58;
k59 = tbl.k59;
k6 = tbl.k6;
k60 = tbl.k60;
k61 = tbl.k61;
k62 = tbl.k62;
k63 = tbl.k63;
k64 = tbl.k64;
k65 = tbl.k65;
k66 = tbl.k66;
k67 = tbl.k67;
k69 = tbl.k69;
k7 = tbl.k7;
k70 = tbl.k70;
k71 = tbl.k71;
k72 = tbl.k72;
k73 = tbl.k73;
k74 = tbl.k74;
k75 = tbl.k75;
k76 = tbl.k76;
k77 = tbl.k77;
k78 = tbl.k78;
k79 = tbl.k79;
k8 = tbl.k8;
k9 = tbl.k9;
k_IFN = tbl.k_IFN;
k_IKK = tbl.k_IKK;
k_IKKe_TBK1 = tbl.k_IKKe_TBK1;
k_IRF3_IKKe_TBK1 = tbl.k_IRF3_IKKe_TBK1;
k_MAVS = tbl.k_MAVS;
k_RIGI = tbl.k_RIGI;
k_TFBS_IFNa = tbl.k_TFBS_IFNa;
k_TFBS_IFNb = tbl.k_TFBS_IFNb;
k_TFBS_IFNl = tbl.k_TFBS_IFNl;
k_a = tbl.k_a;
k_act = tbl.k_act;
k_c = tbl.k_c;
k_deph = tbl.k_deph;
k_e = tbl.k_e;
k_en = tbl.k_en;
k_expr_IkBa = tbl.k_expr_IkBa;
k_f = tbl.k_f;
k_inh_p65 = tbl.k_inh_p65;
k_mRNA_IFNb = tbl.k_mRNA_IFNb;
k_mRNA_IFNl = tbl.k_mRNA_IFNl;
k_mRNA_MX1 = tbl.k_mRNA_MX1;
k_mx1_mRNA = tbl.k_mx1_mRNA;
k_r = tbl.k_r;
k_rigi_synt = tbl.k_rigi_synt;
k_t = tbl.k_t;
k_transISG = tbl.k_transISG;
k_transISGn = tbl.k_transISGn;
k_trans_IFNl = tbl.k_trans_IFNl;
k_transp_NFkB = tbl.k_transp_NFkB;
kinhBySOCS = tbl.kinhBySOCS;
ks = tbl.ks;
mu_IFN = tbl.mu_IFN;
mu_IFNl = tbl.mu_IFNl;
mu_IkBa = tbl.mu_IkBa;
mu_V_E = tbl.mu_V_E;
mu_V_I = tbl.mu_V_I;
mu_mRNA_IFNa = tbl.mu_mRNA_IFNa;
mu_mRNA_IFNb = tbl.mu_mRNA_IFNb;
mu_mRNA_IFNl = tbl.mu_mRNA_IFNl;
mu_mRNA_mx1 = tbl.mu_mRNA_mx1;
mu_mx1 = tbl.mu_mx1;
mu_p = tbl.mu_p;
mu_r = tbl.mu_r;
mu_rigi = tbl.mu_rigi;
nSP = tbl.nSP;
rcsat = tbl.rcsat;
tau = tbl.tau;
tau_5 = tbl.tau_5;

%% Clear temporary variables
clear opts tbl

k_en = (0.43/60);          % Rate of virus entry
k_f = (0.031/60);          % Rate of fusion
k_a = (3.6e-08) * 40;      % Effective virus generation
k_e = (6.6e-02)/60;        % Export rate of RCM
k_r = (3.6/60);            % RNA synthesis rate
k_c = (2.6e-03) * 40;      % Formation rate of RCCM
k_t = (23.7/60);           % Protein production rate
tau = (5.8*60);            % Delay constant
a_RC = (0.028/60);         % Initial RC rate
rcsat = 88/2400;                % Carrying capacity (NC) for HCV 

%%
% Define parameter order
    param = struct();

param.Vc2n = Vc2n;
param.Vn2c = Vn2c;
param.b_IRF3 = b_IRF3;
param.b_MAVS = b_MAVS;
param.b_RIGI = b_RIGI;
param.b_prot = b_prot;

param.k1 = k1;
param.k10 = k10;
param.k11 = k11;
param.k12 = k12;
param.k13 = k13;
param.k14 = k14;
param.k15 = k15;
param.k16 = k16;
param.k17 = k17;
param.k18 = k18;
param.k19 = k19;
param.k2 = k2;
param.k20 = k20;
param.k21 = k21;
param.k22 = k22;
param.k23 = k23;
param.k27 = k27;
param.k29 = k29;
param.k3 = k3;
param.k31 = k31;
param.k32 = k32;
param.k34 = k34;
param.k35 = k35;
param.k36 = k36;
param.k37 = k37;
param.k38 = k38;
param.k39 = k39;
param.k4 = k4;
param.k40 = k40;
param.k41 = k41;
param.k42 = k42;
param.k43 = k43;
param.k44 = k44;
param.k45 = k45;
param.k46 = k46;
param.k47 = k47;
param.k48 = k48;
param.k49 = k49;
param.k5 = k5;
param.k50 = k50;
param.k51 = k51;
param.k52 = k52;
param.k53 = k53;
param.k54 = k54;
param.k56 = k56;
param.k57 = k57;
param.k58 = k58;
param.k59 = k59;
param.k6 = k6;
param.k60 = k60;
param.k61 = k61;
param.k62 = k62;
param.k63 = k63;
param.k64 = k64;
param.k65 = k65;
param.k66 = k66;
param.k67 = k67;
param.k69 = k69;
param.k7 = k7;
param.k70 = k70;
param.k71 = k71;
param.k72 = k72;
param.k73 = k73;
param.k74 = k74;
param.k75 = k75;
param.k76 = k76;
param.k77 = k77;
param.k78 = k78;
param.k79 = k79;
param.k8 = k8;
param.k9 = k9;

param.k_IFN = k_IFN;
param.k_IKK = k_IKK;
param.k_IKKe_TBK1 = k_IKKe_TBK1;
param.k_IRF3_IKKe_TBK1 = k_IRF3_IKKe_TBK1;
param.k_MAVS = k_MAVS;
param.k_RIGI = k_RIGI;

param.k_TFBS_IFNa = k_TFBS_IFNa;
param.k_TFBS_IFNb = k_TFBS_IFNb;
param.k_TFBS_IFNl = k_TFBS_IFNl;

param.k_a = k_a;
param.k_act = k_act;
param.k_c = k_c;
param.k_deph = k_deph;
param.k_e = k_e;
param.k_en = k_en;
param.k_expr_IkBa = k_expr_IkBa;
param.k_f = k_f;
param.k_inh_p65 = k_inh_p65;

param.k_mRNA_IFNb = k_mRNA_IFNb;
param.k_mRNA_IFNl = k_mRNA_IFNl;
param.k_mRNA_MX1 = k_mRNA_MX1;
param.k_mx1_mRNA = k_mx1_mRNA;

param.k_r = k_r;
param.k_rigi_synt = k_rigi_synt;
param.k_t = k_t;

param.k_transISG = k_transISG;
param.k_transISGn = k_transISGn;
param.k_trans_IFNl = k_trans_IFNl;
param.k_transp_NFkB = k_transp_NFkB;

param.kinhBySOCS = kinhBySOCS;
param.ks = ks;

param.mu_IFN = mu_IFN;
param.mu_IFNl = mu_IFNl;
param.mu_IkBa = mu_IkBa;
param.mu_V_E = mu_V_E;
param.mu_V_I = mu_V_I;

param.mu_mRNA_IFNa = mu_mRNA_IFNa;
param.mu_mRNA_IFNb = mu_mRNA_IFNb;
param.mu_mRNA_IFNl = mu_mRNA_IFNl;
param.mu_ISG_RNA = mu_mRNA_mx1;

param.mu_mx1 = mu_mx1;
param.mu_p = mu_p;
param.mu_r = mu_r;
param.mu_rigi = mu_rigi;

param.nSP = nSP;
param.rcsat = rcsat;
param.tau = tau;
param.tau_5 = tau_5;

param.a_RC = a_RC;
param.deg = deg;
param.degARCBySOCS = degARCBySOCS;
param.degRecBySOCS = degRecBySOCS;

param.M = 500;
param.Omega = 50;
param.gamma_RIGI = 1;
param.threshold = 500;
param.n = 2;

save('param_HCV.mat', 'param')


