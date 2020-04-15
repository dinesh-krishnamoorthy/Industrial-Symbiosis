function [sys,F,par] = CentralizedSimulator(par)

import casadi.*

global pO pGL ts

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
rho_ro = sum(rho_o)/2;

mu_oil = 1*0.001; % 1cP oil viscosity

A_w = pi.*(D_w/2).^2;
A_bh = pi.*(D_bh/2).^2;
V_a = L_a.*(pi.*(D_a/2).^2 - pi.*(D_w/2).^2);

% differential states                   Index in solution (sol.x)
m_ga = MX.sym('m_ga',n_w); % 1-6        % 1-6
m_gt = MX.sym('m_gt',n_w); % 7-12       % 7-12
m_ot = MX.sym('m_ot',n_w); % 13-18      % 13-18

% Algebraic states
p_ai = MX.sym('p_ai',n_w);      % 1-6   % 19-24
p_wh = MX.sym('p_wh',n_w);      % 7-12  % 25-30
p_wi = MX.sym('p_wi',n_w);      % 13-18 % 31-36
p_bh = MX.sym('p_bh',n_w);      % 19-24 % 37-42
rho_ai = MX.sym('rho_ai',n_w);  % 25-30 % 43-48
rho_m = MX.sym('rho_m',n_w);    % 31-36 % 49-54
w_iv = MX.sym('w_iv',n_w);      % 37-42 % 55-60
w_pc = MX.sym('w_pc',n_w);      % 43-48 % 61-66
w_pg = MX.sym('w_pg',n_w);      % 49-54 % 67-72
w_po = MX.sym('w_po',n_w);      % 55-60 % 73-78
w_ro = MX.sym('w_ro',n_w);      % 61-66 % 79-84
w_rg = MX.sym('w_rg',n_w);      % 67-72 % 85-90

% control input
w_gl = MX.sym('w_gl',n_w);              % 91-96 

% parameters
p_res = MX.sym('p_res',n_w);
PI = MX.sym('PI',n_w);
GOR = MX.sym('GOR',n_w);      % the only time varying parameter
p_m = MX.sym('p_m',1);
T_a = MX.sym('T_a',n_w);
T_w = MX.sym('T_w',n_w);
R = par.R;
Mw = par.Mw;

f1 = -p_ai.*1e5 + ((R.*T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3) + (Mw./(R.*T_a).*((R.*T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3)).*9.81.*H_a;
f2 = -p_wh.*1e5 + ((R.*T_w./Mw).*(m_gt.*1e3./(L_w.*A_w + L_bh.*A_bh - m_ot.*1e3./rho_o))) - ((m_gt.*1e3+m_ot.*1e3 )./(L_w.*A_w)).*9.81.*H_w/2;
f3 = -p_wi.*1e5 + (p_wh.*1e5 + 9.81./(A_w.*L_w).*max(0,(m_ot.*1e3+m_gt.*1e3-rho_o.*L_bh.*A_bh)).*H_w + 128.*mu_oil.*L_w.*w_pc./(3.14.*D_w.^4.*((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*T_w.*m_gt.*1e3)));
f4 = -p_bh.*1e5 + (p_wi.*1e5 + rho_o.*9.81.*H_bh + 128.*mu_oil.*L_bh.*w_ro./(3.14.*D_bh.^4.*rho_o));
f5 = -rho_ai.*1e2 +(Mw./(R.*T_a).*p_ai.*1e5);
f6 = -rho_m.*1e2 + ((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*T_w.*m_gt.*1e3);
f7 = -w_iv + C_iv.*sqrt(rho_ai.*1e2.*(p_ai.*1e5 - p_wi.*1e5));
f8 = -w_pc + 1.*C_pc.*sqrt(rho_m.*1e2.*(p_wh.*1e5 - p_m.*1e5));
f9 = -w_pg + (m_gt.*1e3./max(1e-3,(m_gt.*1e3+m_ot.*1e3))).*w_pc;
f10 = -w_po + (m_ot.*1e3./max(1e-3,(m_gt.*1e3+m_ot.*1e3))).*w_pc;
f11 = -w_ro + PI.*1e-6.*(p_res.*1e5 - p_bh.*1e5);
f12 = -w_rg.*1e-1 + GOR.*w_ro;

df1 = (w_gl - w_iv).*1e-3; 
df2 = (w_iv + w_rg.*1e-1 - w_pg).*1e-3; 
df3 = (w_ro - w_po).*1e-3; 

% Form the DAE system
diff = vertcat(df1,df2,df3);
alg = vertcat(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12);
x_var = vertcat(m_ga,m_gt,m_ot);
z_var = vertcat(p_ai,p_wh,p_wi,p_bh,rho_ai,rho_m,w_iv,w_pc,w_pg,w_po,w_ro,w_rg);
p_var = vertcat(w_gl,GOR,p_m);
u_var = w_gl;
d_var = vertcat(GOR,p_m);
y_var = vertcat(p_wh, p_bh, w_po, w_pg, w_ro, w_rg);

L = -pO*sum(w_po) + pGL*sum(w_gl); 

alg = substitute(alg,p_res,par.p_res);
alg = substitute(alg,PI,par.PI);
alg = substitute(alg,T_a,par.T_a);
alg = substitute(alg,T_w,par.T_w);

dae = struct('x',x_var,'z',z_var,'p',p_var,'ode',diff,'alg',alg,'quad',L); ...
    % der(m_tot) = w_in - w_out;
opts = struct('tf',ts);

% create IDAS integrator
F = integrator('F','idas',dae,opts);

%%

sys.x = x_var;
sys.z = z_var;
sys.u = p_var;
sys.d = d_var;
sys.diff = diff;
sys.alg = alg;
sys.L = L;
sys.y = y_var;
sys.f = [];

sys.nlcon = [];
sys.lb =  [];
sys.ub = [];
