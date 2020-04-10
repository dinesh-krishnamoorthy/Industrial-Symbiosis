function [F,x_init,par] = PlantSimulatorModel(Ts,par1,par2)

% Geneerate the plant simulator using well parameters of the different
% companies. 

par = cell2struct(cellfun(@vertcat,struct2cell(par1),struct2cell(par2),...
    'uni',0),fieldnames(par1),1);

par.n_w = sum(par.n_w); % Total number of wells

% Create a centralized plant simulator. 
[sys,F,~,par] = GasLiftModel(par,Ts);

% Return steady-state initial state for the following initial conditions.
u_init = [1,1,1,1,1,1]'; 
d_init = vertcat(par.GOR,20);
x_init = solveODE(sys,par,d_init,u_init);

par.ind1 = [1:3,7:9,13:15];
par.ind2 = [4:6,10:12,16:18];
