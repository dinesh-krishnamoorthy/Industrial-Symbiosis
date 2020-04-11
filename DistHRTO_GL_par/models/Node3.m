function [sys,par] = Node3(par)

import casadi.*
global pGL

if nargin<1
    % Default values
    par.w_pg_max  = 27.5;
    par.w_in_max  = 8.5;
end

w_gl_FPSO = MX.sym('w_gl_FPSO',1);  %value is -ve to represent production
w_pg_FPSO = MX.sym('w_po_FPSO',1);  %value is -ve to represent consump

sys.x = [w_gl_FPSO;w_pg_FPSO];
sys.u = [];
sys.diff = [];
sys.d = [];
sys.L = 0;

sys.nlcon = [];
sys.lb =  [];
sys.ub = [];

par.dx0 = [-8;-27];
par.lbx = [-par.w_in_max;-par.w_pg_max];
par.ubx = [0;0];
par.u0 = [];
par.lbu = [];
par.ubu = [];

