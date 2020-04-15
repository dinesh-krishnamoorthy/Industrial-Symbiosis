function [Ju_hat,Gu_hat] = GradientEstimator(EKF,x_hat,u,d)
% Written by Krishnamoorthy, Oct 2018.
% email - dineshk@ntnu.no

uEXO = vertcat(u,d);

A = full(EKF.JacAx(x_hat,uEXO));
B = full(EKF.JacBu(x_hat,uEXO));
C = full(EKF.JacJx(x_hat,uEXO));
D = full(EKF.JacJu(x_hat,uEXO));
Ju_hat = (-C*(A\B) + D)';

C2 = full(EKF.JacGx(x_hat,uEXO));
D2 = full(EKF.JacGu(x_hat,uEXO));
Gu_hat = (-C2*(A\B) + D2)';

end