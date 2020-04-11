function par = WellParameters_Company2

par.n_w = 3;

par.L_w = [1500;1500;1500];
par.H_w = [1000;1000;1000];
par.D_w = [0.121;0.121;0.121];

par.L_bh = [500;500;500];
par.H_bh = [500;500;500];
par.D_bh = [0.121;0.121;0.121];

par.L_a = par.L_w;
par.H_a = par.H_w;
par.D_a = [0.189;0.189;0.189];

par.mu_oil = 1*0.001.*ones(3,1); %1cP
par.rho_o = [8;8.2;8.05].*1e2;
par.C_iv = [0.1e-3;0.1e-3;0.1e-3];
par.C_pc = [2e-3;2e-3;2e-3];

par.GOR = [0.108;0.115;0.102];
par.p_res = [160;155;155];
par.PI = [7;7;7].*0.5;
par.T_a = [28+273;28+273;28+273];
par.T_w = [32+273;32+273;32+273];

par.GOR_var = [0.02;0.02;0.02];

par.tf = 1;
par.tSim = 1;
par.nIter = 60;
