function [data,init] = PID_control(data,init)

global sim_k

% --------- Well 1 ---------
tauC = 200;
tauI = min(4*tauC,800);
Kp1 = 80/(0.439*tauC);
Ki1 = Kp1/tauI;

init.err(1)= -data.AL(1,sim_k);

data.u(1,sim_k+1) = max(0.1,min(12,  data.u(1,sim_k) + ...
    (Kp1*init.err(1) + Ki1*init.err(1) - Kp1*init.err0(1))));
init.err0(1) = init.err(1);

% --------- Well 2 ---------
tauC = 200;
tauI = min(4*tauC,800);
Kp1 = 80/(0.29*tauC);
Ki1 = Kp1/tauI;

init.err(2) =-data.AL(2,sim_k);

data.u(2,sim_k+1) = max(0.1,min(12, data.u(2,sim_k) +...
    (Kp1*init.err(2) + Ki1*init.err(2) - Kp1*init.err0(2))));
init.err0(2) = init.err(2);

% --------- Well 3 ---------
tauC = 200;
tauI = min(4*tauC,800);
Kp1 = 80/(0.439*tauC);
Ki1 = Kp1/tauI;

init.err(3)= -data.AL(3,sim_k);

data.u(3,sim_k+1) = max(0.1,min(12, data.u(3,sim_k) +...
    (Kp1*init.err(3) + Ki1*init.err(3) - Kp1*init.err0(3))));
init.err0(3) = init.err(3);