function par = Param_6wells(Ind)

par.Ind = Ind;
par.n_w = length(Ind);

par.L_w = [1500;1500;1500;1500;1500;1500]; 
par.H_w = [1000;1000;1000;1000;1000;1000]; 
par.D_w = [0.121;0.121;0.121;0.121;0.121;0.121];
 
par.L_bh = [500;500;500;500;500;500]; 
par.H_bh = [500;500;500;500;500;500]; 
par.D_bh = [0.121;0.121;0.121;0.121;0.121;0.121]; 
 
par.L_a = par.L_w; 
par.H_a = par.H_w; 
par.D_a = [0.189;0.189;0.189;0.189;0.189;0.189]; 
 
par.rho_o = [8;8;7.9;8;8.2;8.05].*1e2; 
par.C_iv = [0.1e-3;0.1e-3;0.1e-3;0.1e-3;0.1e-3;0.1e-3]; 
par.C_pc = [2e-3;2e-3;2e-3;2e-3;2e-3;2e-3]; 
 
par.GOR = [0.1;0.12;0.09;0.108;0.115;0.102]; 
par.p_m = [20;20;20;20;20;20]; 
par.p_res = [150;155;155;160;155;155]; 
par.PI = [7;7;7;7;7;7].*0.5; 
par.T_a = [28+273;28+273;28+273;28+273;28+273;28+273]; 
par.T_w = [32+273;32+273;32+273;32+273;32+273;32+273]; 
 
par.GOR_var = [0.02;0.02;0.02;0.02;0.02;0.02]; 
  
par.R = 8.314; 
par.Mw = 20e-3; 


%%
par.L_w = par.L_w(Ind);
par.H_w = par.H_w(Ind);
par.D_w = par.D_w(Ind);
par.L_bh = par.L_bh(Ind);
par.H_bh = par.H_bh(Ind);
par.D_bh = par.D_bh(Ind);
par.L_a = par.L_a(Ind);
par.H_a = par.H_a(Ind);
par.rho_o = par.rho_o(Ind);
par.C_iv = par.C_iv(Ind);
par.D_a = par.D_a(Ind);
par.C_pc = par.C_pc(Ind);
par.GOR = par.GOR(Ind);
par.p_m = par.p_m(Ind);
par.p_res = par.p_res(Ind);
par.PI = par.PI(Ind);
par.T_a = par.T_a(Ind);
par.T_w = par.T_w(Ind);
par.GOR_var = par.GOR_var(Ind);