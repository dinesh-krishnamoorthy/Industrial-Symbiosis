function EKF = EKF_setup(sys,P,Q,R)
import casadi.*

xEKF = vertcat(sys.x,sys.d);

EKF.f = Function('f_EKF',{xEKF,sys.u},{sys.x_k1},{'x','p'},{'xdot'});
EKF.JacFx = Function('JacFx',{xEKF,sys.u},{jacobian(sys.x_k1,xEKF)});

EKF.h = Function('h_EKF',{xEKF,sys.u},{sys.y});
EKF.JacHx = Function('JacHx',{xEKF,sys.u},{jacobian(sys.y,xEKF)});

p_var = vertcat(sys.u,sys.d);
EKF.JacAx = Function('JacAx',{sys.x,p_var},{jacobian(sys.diff,sys.x)});
EKF.JacBu = Function('JacBu',{sys.x,p_var},{jacobian(sys.diff,sys.u)});
EKF.JacJx = Function('JacJx',{sys.x,p_var},{jacobian(sys.L,sys.x)});
EKF.JacJu = Function('JacJu',{sys.x,p_var},{jacobian(sys.L,sys.u)});

EKF.JacGx = Function('JacGx',{sys.x,p_var},{jacobian(sys.g,sys.x)});
EKF.JacGu = Function('JacGu',{sys.x,p_var},{jacobian(sys.g,sys.u)});

if nargin<2
    EKF.Pk = 1e3.*eye(numel(xEKF));
    EKF.Qk = 1e3.*eye(numel(xEKF));
    EKF.Rk = 1e0.*eye(numel(sys.y));
else
    EKF.Pk = P;
    EKF.Qk = Q;
    EKF.Rk = R;
end




