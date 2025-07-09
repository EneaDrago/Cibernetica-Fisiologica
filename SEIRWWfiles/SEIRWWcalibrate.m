function [params,J] = SEIRWWcalibrate(YC,YW,C,params)
% Parameter estimation

%Rate E -> I
params.alpha = .4433;

%Initial rate S -> I
params.beta = .44;

%Rate I to R (tau1 in SEIR-ICU model)
params.tau = .32; 

%State noise coefficient (model error)
params.modelErrorC = 4^2;

%Initial error variance of beta
params.S_beta = .15^2;

%Variance of daily change of beta (initially)
params.Q_beta0 = .05^2;

%After 1st month
params.Q_beta1 = .005^2;

%Estimate the initial sizes of E and I compartments
params.E_init = params.darkNumber(1,1)/params.alpha * (1 + mean(YC(1:5)));
params.I_init = params.darkNumber(1,1)/params.tau * (1 + mean(YC(1:5)));

params.varE_init = (params.E_init/2)^2;
params.varI_init = (params.I_init/2)^2;

%Find initial point by a simpler optimisation
params.gamma = 2;
params.WWexp = .7;
cost = @(x)(paramFit(params,YC,YW,C,x(1),x(2),-1));
Acon = [1 0; -1 0; 0 1; 0 -1];
bcon = [4 -.2 1 -.4];
xopt = fmincon(cost,[params.gamma params.WWexp],Acon,bcon);
[~, nuInit, ~] = cost(xopt);
params.gamma = xopt(1);
params.WWexp = xopt(2);

disp('Initial point found')


%Estimate gamma, nu, and the exponent in WW-transformation
cost = @(x)(paramFit(params,YC,YW,C,x(1),x(2),x(3)));
Acon = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 -1];
bcon = [4 -.2 1 -.4 0];
xopt = fmincon(cost,[params.gamma params.WWexp nuInit],Acon,bcon);
[J, ~, RW0] = cost(xopt);
params.gamma = xopt(1);
params.WWexp = xopt(2);
params.nu = xopt(3);
params.RW0 = RW0;

disp('Parameter estimation complete')