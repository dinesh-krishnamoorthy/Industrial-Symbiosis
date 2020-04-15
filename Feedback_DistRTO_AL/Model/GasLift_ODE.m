function [sys,init] = GasLift_ODE(par)

import casadi.*

global pO pGL tSim
n_w = par.n_w;

L_w = par.L_w;
H_w = par.H_w;
D_w = par.D_w;

L_bh = par.L_bh;
H_bh = par.H_bh;
D_bh = par.D_bh;

L_a = par.L_a;
H_a = par.H_a;
D_a = par.D_a;

rho_o = par.rho_o;
C_iv = par.C_iv;
C_pc = par.C_pc;  

mu_oil = 1*0.001; % 1cP oil viscosity

A_w = pi.*(D_w/2).^2;
A_bh = pi.*(D_bh/2).^2;
V_a = L_a.*(pi.*(D_a/2).^2 - pi.*(D_w/2).^2);

Mw = par.Mw;
R = par.R;

% differential states
m_ga = MX.sym('m_ga',n_w); % 1-2
m_gt = MX.sym('m_gt',n_w); % 3-4
m_ot = MX.sym('m_ot',n_w); % 5-6

% control input
w_gl = MX.sym('w_gl',n_w);
GOR = MX.sym('GOR',n_w);

% algebraic equations used for substitution in the ODE model
p_m = par.p_m;
p_ai = 1e-5.*(((R.*par.T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3) + (Mw./(R.*par.T_a).*((R.*par.T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3)).*9.81.*H_a);
p_wh = 1e-5.*(((R.*par.T_w./Mw).*(m_gt.*1e3./(L_w.*A_w + L_bh.*A_bh - m_ot.*1e3./rho_o))) - ((m_gt.*1e3+m_ot.*1e3 )./(L_w.*A_w)).*9.81.*H_w/2);
rho_ai = 1e-2.*(Mw./(R.*par.T_a).*p_ai.*1e5);
rho_m = 1e-2.*(((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*par.T_w.*m_gt.*1e3)); 
w_pc = C_pc.*sqrt(rho_m.*1e2.*(p_wh.*1e5 - p_m.*1e5));
w_pg = (m_gt.*1e3./(m_gt.*1e3+m_ot.*1e3)).*w_pc;
w_po = (m_ot.*1e3./(m_gt.*1e3+m_ot.*1e3)).*w_pc;
p_wi = 1e-5.*((p_wh.*1e5 + 9.81./(A_w.*L_w).*max(0,(m_ot.*1e3+m_gt.*1e3-rho_o.*L_bh.*A_bh)).*H_w + 128.*mu_oil.*L_w.*w_pc./(3.14.*D_w.^4.*((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*par.T_w.*m_gt.*1e3))));
p_bh = 1e-5.*(p_wi.*1e5 + rho_o.*9.81.*H_bh + 128.*mu_oil.*L_bh.*w_po./(3.14.*D_bh.^4.*rho_o));
w_iv = C_iv.*sqrt(rho_ai.*1e2.*(p_ai.*1e5 - p_wi.*1e5));
w_ro = (par.PI).*1e-6.*(par.p_res.*1e5 - p_bh.*1e5);% w_ro = (-IPR.a + sqrt(IPR.a.^2+4*IPR.b.*(p_res-p_bh).*1e5))./(2.*IPR.b);
w_rg = 1e1.*GOR.*w_ro;


% differential equations
df1C = (w_gl - w_iv).*1e-3;
df2C = (w_iv + w_rg.*1e-1 - w_pg).*1e-3;
df3C = (w_ro - w_po).*1e-3;

% Concatenate the differential and algebraic equations
diff = vertcat(df1C,df2C,df3C);

% Difference equations - Descretized model for EKF
df1 = m_ga + tSim.*(w_gl - w_iv).*1e-3;
df2 = m_gt + tSim.*(w_iv + w_rg.*1e-1 - w_pg).*1e-3;
df3 = m_ot + tSim.*(w_ro - w_po).*1e-3;
df4 = GOR ;

% Concatenate the differential and algebraic equations
diff_EKF = vertcat(df1,df2,df3,df4);

% concatenate the differential and algebraic states
x_var = vertcat(m_ga,m_gt,m_ot);
u_var = vertcat(w_gl);
d_var = vertcat(GOR);
x_EKF = vertcat(m_ga,m_gt,m_ot,GOR);
y_EKF = vertcat(p_wh, p_bh, w_po, w_pg, w_ro, w_rg);

z_var = vertcat(p_ai,p_wh,p_wi,p_bh,rho_ai,rho_m,w_iv,w_pc,w_pg,w_po,w_ro,w_rg);
alg = Function('alg',{x_var,u_var,d_var},{z_var});

L = -pO*sum(w_po) + pGL*sum(w_gl); 

sys.x = x_var;
sys.u = u_var;
sys.d = d_var;
sys.x_k1 = diff_EKF;
sys.diff = diff;
sys.L = L;
sys.y = y_EKF;
sys.alg = alg;

sys.nlcon = [];
sys.lb =  [];
sys.ub = [];

p_wh_Index = n_w+1:2*n_w;
p_bh_Index = 3*n_w+1:4*n_w;
w_po_Index = 9*n_w+1:10*n_w;
w_pg_Index = 8*n_w+1:9*n_w;
w_ro_Index = 10*n_w+1:11*n_w;
w_rg_Index = 11*n_w+1:12*n_w;
        
sys.yIndex = [p_wh_Index, p_bh_Index,w_po_Index,w_pg_Index,w_ro_Index, w_rg_Index];

%% 

par.dx0 = [1.32.*ones(n_w,1);0.8.*ones(n_w,1);6.*ones(n_w,1)];
par.lbx = [0.01.*ones(n_w,1);0.01.*ones(n_w,1);0.01.*ones(n_w,1)];
par.ubx = [1e7.*ones(n_w,1);1e7.*ones(n_w,1);1e7.*ones(n_w,1)];
par.u0 = 1.*ones(n_w,1);
par.lbu = 0.9.*ones(n_w,1);
par.ubu = 10.*ones(n_w,1);

init.d = par.GOR;
init.u = ones(n_w,1);
[init.xf,init.zf,init.exitflag] = solveODE(sys,par,init.d,init.u);
