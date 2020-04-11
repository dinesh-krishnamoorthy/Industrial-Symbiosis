function [MPCsolver,MPC_par] = SetpointNMPC(sys,par)

import casadi.*

% Direct Collocation - Polynomials
d = 3;
[B,C,D] = DirectCollocation(d);

%% Build NLP solver

nx = numel(sys.x);
nu = numel(sys.u);
nd = numel(sys.d);

assert(numel(sys.x)==numel(par.lbx),'Dimension mismatch in x.')
assert(numel(sys.x)==numel(par.ubx),'Dimension mismatch in x.')
assert(numel(sys.x)==numel(par.dx0),'Dimension mismatch in x.')
assert(numel(sys.u)==numel(par.lbu),'Dimension mismatch in u.')
assert(numel(sys.u)==numel(par.ubu),'Dimension mismatch in u.')
assert(numel(sys.u)==numel(par.u0),'Dimension mismatch in u.')

% empty nlp
w = {};     
w0 = [];        
lbw = [];         
ubw = [];

g = {};     
lbg = [];       
ubg = [];

J = 0;

% initial conditions for each scenario
X0 = MX.sym('X0',nx);  
w = {w{:}, X0};         
lbw = [lbw; par.lbx];         
ubw = [ubw; par.ubx];
w0 = [w0; par.dx0];

X0_par = MX.sym('X0_par',nx);

% initial conditions
g = {g{:},X0 - X0_par};  
lbg = [lbg;zeros(nx,1)];  
ubg = [ubg;zeros(nx,1)];

U0 = MX.sym('U0',nu);

d_est = MX.sym('GOR_est',par.n_w+1);
w_pg_SP = MX.sym('w_pg_SP',par.n_w);

% Formulate NLP
Xk = X0;        
Uk_prev = U0;

for k = 0:par.Horizon_Samples-1
    
    Uk = MX.sym(['U_' num2str(k)],nu);
    
    Upar = vertcat(Uk,d_est);
    
    w = {w{:},Uk};
    lbw = [lbw;par.lbu];
    ubw = [ubw;par.ubu];
    w0 = [w0;par.u0];
    
    Xkj = {};
    Zkj = {};
    
    % New Collocation variables that are handles for solver
    for j = 1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)],nx);
       
        w = {w{:},Xkj{j}};
        lbw = [lbw;par.lbx];
        ubw = [ubw;par.ubx];
        w0 = [w0;par.dx0];
    end
    
    % Loop over collocation points
    Xk_end = D(1)*Xk;
    
    for j = 1:d
        % Expression for the state derivative at the collocation point
        xp = C(1,j+1)*Xk;  % helper state
        for r = 1:d
            xp = xp + C(r+1,j+1)*Xkj{r};
        end
        [fj] =  sys.f(Xkj{j},Upar);
        
        g = {g{:},par.tf*fj-xp};  % dynamics and algebraic constraints
        lbg = [lbg;zeros(nx,1)];
        ubg = [ubg;zeros(nx,1)];
        
        % Add contribution to the end states
        Xk_end = Xk_end + D(j+1)*Xkj{j};
    end
    
    % Slack variables for gas capacity constraints
    s{j} = MX.sym(['s_' num2str(k)],1);
    
    w = {w{:},s{j}};
    lbw = [lbw;0];
    ubw = [ubw;1];
    w0 = [w0;0];
    
    
    Ykj = sys.y_model(Xkj{j},Upar);
    
    J = J + (Ykj(3*par.n_w+1:4*par.n_w)-w_pg_SP)'*(Ykj(3*par.n_w+1:4*par.n_w)-w_pg_SP) +...
        40*s{j} +  1*sum((Uk_prev - Uk).^2);
    
    Uk_prev = MX.sym(['Uprev_' num2str(k+1)],nu);
    Uk_prev = Uk;
    
    % Max produced gas constraint
    g = {g{:},sum(Ykj(3*par.n_w+1:4*par.n_w))-s{j}};
    lbg = [lbg;0];
    ubg = [ubg;27.5];
    
    % Max gas lift constraint
    g = {g{:},sum(Uk)};
    lbg = [lbg;0];
    ubg = [ubg;9.5];
    
    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], nx);
    w = {w{:},Xk};
    lbw = [lbw;par.lbx];
    ubw = [ubw;par.ubx];
    w0 = [w0; par.dx0];
    
    % Shooting Gap constraint
    g = {g{:},Xk_end-Xk};
    lbg = [lbg;zeros(nx,1)];
    ubg = [ubg;zeros(nx,1)];
    
end

%% create NLP solver

opts = struct('warn_initial_bounds',false,'print_time',false, ...
    'ipopt',struct('print_level',1));

nlp = struct('x',vertcat(w{:}),'p',vertcat(w_pg_SP,d_est,X0_par,U0),...
    'f',J,'g',vertcat(g{:}));

MPCsolver = nlpsol('solver','ipopt',nlp,opts);

MPC_par.w0 = w0;
MPC_par.lbw = lbw;
MPC_par.ubw = ubw;
MPC_par.lbg = lbg;
MPC_par.ubg = ubg;
MPC_par.d = d; % no. of collocation points
