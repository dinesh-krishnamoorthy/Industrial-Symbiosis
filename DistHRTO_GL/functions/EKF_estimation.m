function [xk_hat,Pk] = EKF_estimation(EKF,ymeas,uEKF,xk_hat,Pk)

nxEKF = numel(xk_hat);

Fj = full(EKF.JacFx(xk_hat,uEKF));

xk_hat_1 =  full(EKF.f(xk_hat,uEKF));
Pk_1 = Fj*Pk*Fj' + EKF.Qk;

Hj = full(EKF.JacHx(xk_hat_1,uEKF));
ek = full(ymeas - EKF.h(xk_hat_1,uEKF));
Sk = Hj*Pk_1*Hj' + EKF.Rk;
Kk = (Pk_1*Hj')/(Sk);
xk_hat = xk_hat_1 + Kk*ek;
Pk = (eye(nxEKF) - Kk*Hj)*Pk_1;


end