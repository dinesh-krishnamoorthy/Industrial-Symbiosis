function [solver,par] = AugmentedLagrangian_sharing(sys,par,opts)

% Function that computes the steady-state optimum
% Written by Dinesh Krishnamoorthy, Jul 2019, NTNU

import casadi.*

if nargin<8
    opts = struct('warn_initial_bounds',false, ...
        'print_time',false, ...
        'ipopt',struct('print_level',1)...
        );
end

assert(numel(sys.x)==numel(par.lbx),'Dimension mismatch in x.')
assert(numel(sys.x)==numel(par.ubx),'Dimension mismatch in x.')
assert(numel(sys.x)==numel(par.dx0),'Dimension mismatch in x.')
assert(numel(sys.u)==numel(par.lbu),'Dimension mismatch in u.')
assert(numel(sys.u)==numel(par.ubu),'Dimension mismatch in u.')
assert(numel(sys.u)==numel(par.u0),'Dimension mismatch in u.')
assert(~isempty(sys.g),'No Coupling constraints found. Use SSopt instead.')

lambda = MX.sym('lambda',numel(sys.g));
g_k = MX.sym('g_k',numel(sys.g));
g_avg = MX.sym('g_avg',numel(sys.g));
rho = MX.sym('rho');

w = {};
w0 = [];
lbw = [];
ubw = [];

g = {};
lbg = [];
ubg = [];

w = {w{:},sys.x,sys.u};
lbw = [lbw;par.lbx;par.lbu];
ubw = [ubw;par.ubx;par.ubu];
w0 = [w0;par.dx0;par.u0];

J = sys.L + lambda'*sys.g + rho/2*norm(sys.g -  (g_k - g_avg)).^2; % Economic objective + coupling constraints

if ~isempty(sys.diff)
    g = {g{:},vertcat(sys.diff)};
    lbg = [lbg;zeros(numel(sys.diff),1)];
    ubg = [ubg;zeros(numel(sys.diff),1)];
end

if ~isempty(sys.nlcon)
    assert(numel(sys.nlcon)==numel(sys.lb))
    assert(numel(sys.nlcon)==numel(sys.ub))
    
    g = {g{:},sys.nlcon};
    lbg = [lbg;sys.lb];
    ubg = [ubg;sys.ub];
end

nlp = struct('x',vertcat(w{:}),'p',vertcat(sys.d,lambda,rho,g_k,g_avg),'f',J,'g',vertcat(g{:}));
solver = nlpsol('solver','ipopt',nlp,opts);

par.w0 = w0;
par.lbw = lbw;
par.ubw = ubw;
par.lbg = lbg;
par.ubg = ubg;
par.nlp = nlp;




