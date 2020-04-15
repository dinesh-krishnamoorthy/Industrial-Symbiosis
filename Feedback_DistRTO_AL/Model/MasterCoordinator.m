function public = MasterCoordinator(public)

global sim_k qGLMax rho

public.w_gl3_opt(sim_k) =  max(-qGLMax,...
    -public.lambda_wgl(sim_k)/rho ...
    -public.w_gl1_opt(sim_k)...
    -public.w_gl2_opt(sim_k));

public.residual(sim_k) = ...
      public.w_gl1_opt(sim_k) ...
    + public.w_gl2_opt(sim_k) ...
    + public.w_gl3_opt(sim_k);

public.lambda_wgl(sim_k+1) = ...
    public.lambda_wgl(sim_k) + rho*(public.residual(sim_k));

public.dlambda_wgl(sim_k) = public.lambda_wgl(sim_k+1) - public.lambda_wgl(sim_k);

public.lambda_wpg(sim_k+1) = 0;% lambda_wpg + rho_wpg * ...
%(lambda1_star_Qg + lambda2_star_Qg - par.QgMax);


