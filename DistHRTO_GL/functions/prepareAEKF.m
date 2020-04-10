function EKF = prepareAEKF(sys,Ts)

import casadi.*

nx = numel(sys.x);
nd = numel(sys.d);
ny = length(sys.y);

diff_EKF= [];

for i = 1:nx
    diff_EKF = [diff_EKF; sys.x(i) + Ts*sys.diff(i)];
end
for i = 1:nd
    diff_EKF = [diff_EKF; sys.d(i)];
end

x_EKF = vertcat(sys.x,sys.d);

EKF.f = Function('f_EKF',{x_EKF,sys.u},{diff_EKF},{'x','p'},{'xdot'});
EKF.JacFx = Function('JacFx',{x_EKF,sys.u},{jacobian(diff_EKF,x_EKF)});

EKF.h = Function('h_EKF',{x_EKF,sys.u},{sys.y});
EKF.JacHx = Function('JacHx',{x_EKF,sys.u},{jacobian(sys.y,x_EKF)});

p_var = vertcat(sys.u,sys.d);
EKF.JacAx = Function('JacAx',{sys.x,p_var},{jacobian(sys.diff,sys.x)});
EKF.JacBu = Function('JacBu',{sys.x,p_var},{jacobian(sys.diff,sys.u)});
EKF.JacJx = Function('JacJx',{sys.x,p_var},{jacobian(sys.L, sys.x)});
EKF.JacJu = Function('JacJu',{sys.x,p_var},{jacobian(sys.L, sys.u)});


EKF.Pk = 1e3.*eye(nx+nd);
EKF.Qk = 1e3.*eye(nx+nd);
EKF.Rk = 1e0.*eye(ny);


