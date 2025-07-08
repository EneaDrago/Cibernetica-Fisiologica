function [Yest, err] = SEIR_WW_FWD(X,C,P,initDay,params,maxind)

%X: initial state for projection
%C: daily observation coefficients
%initDay: number of the first day of simulation
%params: the coefficients in the SEIR ODE
%maxInd: length of simulation time

% Set parameters
alpha = params.alpha;
tau = params.tau;
gamma = params.gamma;
nu = params.nu;   
eta = 1;
N = params.N;
CC = params.modelErrorC; 

%Reaction stoichiometry (w.r.t. SEIRreaction-function)
AR = [-1 0 0 0 0 0;
    1 -1 0 0 0 0;
    0 1 -1 0 0 0;
    1 0 0 -1 0 0;
    0 1 0 0 0 0;
    0 0 0 0 1 -1;
    0 0 0 0 0 0];

%Time step = 1/N_step (days)
N_step = 10;

%Number of detected cases today depends linearly on the true number of new
%cases today
Ccase = [0, 0, 0, 0, 1, 0, 0];
Cww = [0, 0, 0, 0, 0, 1, 0];

Yest = zeros(2,maxind);
err = zeros(1,maxind);

% Loop sui giorni da stimare
for jday = 1:maxind
    
    %Reset the "cases today" counter
    X(5) = 0;
    P(5,:) = 0;
    P(:,5) = 0;
    
    %Time loop for one day
    for jt = 1:N_step
        [RR, Jf] = SEIRreaction(X,N,alpha,tau,gamma,nu,eta,1/N_step);
        X = X + AR*RR;
        Q = CC*AR*diag(RR)*AR';
        P = (eye(7) + AR*Jf)*P*(eye(7) + Jf'*AR') + Q;
    end
    
    %Predicted number of daily new cases and wastewater measurement
    Yest(:,jday) = [C(initDay+jday-1)*Ccase; Cww]*X;
    err(jday) = P(5,5).^.5;
end

