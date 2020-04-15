function [data,init] = SubsystemEstimator(EKF,data,public,init)
import casadi.*

global sim_k rho

[init.xk_hat,init.Pk] = EKF_estimation(EKF,data.y(:,sim_k),...
        data.u(:,sim_k),init.xk_hat,init.Pk);
    x_hat = init.xk_hat(1:9);
    d_hat = init.xk_hat(10:12);
    
    [Ju_hat,Gu_hat] = GradientEstimator(EKF,x_hat,data.u(:,sim_k),d_hat);
    
    % Update Ju_hat
    AL_hat = Ju_hat + public.lambda_wgl(sim_k)*Gu_hat(:,1) + ...
        public.lambda_wpg(sim_k)*Gu_hat(:,2) + ...
        rho*(public.residual(sim_k));
    
    data.Ju(:,sim_k) = Ju_hat;
    data.gu1(:,sim_k) = Gu_hat(:,1);
    data.gu2(:,sim_k) = Gu_hat(:,2);
    data.AL(:,sim_k) = AL_hat;
    