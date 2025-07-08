function [R, JR] = SEIRreaction(X,N,alpha,tau1,gamma,nu,eta,dt)

% X1: S(t)
% X2: E(t)
% X3: I(t)
% X4: A(t)
% X5: N(t)
% X6: W(t)
% X7: beta(t)

%Reactions
R = zeros(6,1);
R(1) = dt*X(7)*X(1)*X(3)/N;     % S to E
R(2) = dt*alpha*X(2);           % E to I
R(3) = dt*tau1*X(3);             % I to R
R(4) = dt*gamma*X(4);           % A to void
R(5) = dt*nu*X(4);              % W production
R(6) = dt*eta*X(6);               % W degradation

R = max(R,0);

%Jacobian of R
JR = zeros(6,7);
JR(1,:) = [dt*X(7)*X(3)/N, 0, dt*X(7)*X(1)/N, 0, 0, 0, dt*X(1)*X(3)/N];
JR(2,:) = [0, dt*alpha, 0, 0, 0, 0, 0];
JR(3,:) = [0, 0, dt*tau1, 0, 0, 0, 0];
JR(4,:) = [0, 0, 0, dt*gamma, 0, 0, 0];
JR(5,:) = [0, 0, 0, dt*nu, 0, 0, 0];
JR(6,:) = [0, 0, 0, 0, 0, dt*eta, 0];


%Reaction stoichiometry w.r.t. reactions defined here
% AR = [-1 0 0 0 0 0;
%     1 -1 0 0 0 0;
%     0 1 -1 0 0 0;
%     1 0 0 -1 0 0;
%     0 1 0 0 0 0;
%     0 0 0 0 1 -1;
%     0 0 0 0 0 0];




