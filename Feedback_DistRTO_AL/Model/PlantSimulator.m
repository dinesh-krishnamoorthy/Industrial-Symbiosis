function [data1,data2,public,init] = PlantSimulator(F,init,data1,data2,public,par,sys)

import casadi.*

global sim_k qGLMax QgMax

GOR_real = par.GOR_val_SS;

% ----------- vary GOR conditions in the wells -----------
if sim_k >3*3600
    GOR_real(3) = 0.1100;
end
if sim_k >4*3600
    GOR_real(5) = 0.1320;
end
if sim_k>2*3600 && sim_k<=7*3600
    qGLMax = 6;
end
if sim_k>7*3600 && sim_k<=8*3600
    qGLMax = 4;
end
if sim_k>8*3600
    qGLMax = 7;
end
public.qGLMax(sim_k) = qGLMax;
public.QgMax(sim_k) = QgMax;


u_in = [init.u;GOR_real;20];

meas.GOR_real(:,sim_k) = GOR_real;

% Simulator using IDAS integrator
Fk = F('x0',init.xf,'z0',init.zf,'p',u_in);
init.xf = full(Fk.xf);
init.zf = full(Fk.zf);
J_real(sim_k) = full(Fk.qf);

n_w = par.n_w;
n_w1 = 3;

% Measured Data
ymeas = init.zf(sys.yIndex);
data1.y(:,sim_k+1) = [ymeas(1:n_w1);ymeas(n_w+1:n_w+n_w1);ymeas(2*n_w+1:2*n_w+n_w1);ymeas(3*n_w+1:3*n_w+n_w1);ymeas(25:27);ymeas(31:33)];
data2.y(:,sim_k+1) = [ymeas(n_w1+1:n_w);ymeas(n_w+n_w1+1:2*n_w);ymeas(2*n_w+n_w1+1:3*n_w);ymeas(3*n_w+n_w1+1:4*n_w);ymeas(28:30);ymeas(34:36)];

public.w_gl1_opt(sim_k+1) = sum(init.u(1:3));
public.w_pg1_opt(sim_k+1) = sum(data1.y(10:12,sim_k+1));

public.w_gl2_opt(sim_k+1) = sum(init.u(4:6));
public.w_pg2_opt(sim_k+1) = sum(data2.y(10:12,sim_k+1));
