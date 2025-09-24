%% KYA314 Assignment 3, Q3c

clear;
clc;
format long;

%% Define parameters

alpha = 0.5;
beta = 0.2;
delta = 0.2;
K = 1;

%% Define the Lotka-Volterra equations

lotkaVolterra = @(x)([alpha * x(1) * (1 - x(1) / K) - beta * x(1) * x(2); beta * x(1) * x(2) - delta * x(2)]);

%% Define the Jacobian function

lvJacobian = @(x)([alpha - 2 * alpha * x(1) / K - beta * x(2), -1 * beta * x(1); beta * x(2), beta * x(1) - delta]);

%% Run MySolve to find the equilibria

[x1, conv1, J1] = MySolve(lotkaVolterra, [1.1;1.1], lvJacobian, 10^-6, 1000);

x1

% Thus [x,y]=[1,0] is the prey-only equilibrium and the coexistence
% equilibrium, as with these values delta/beta=1, and alpha/beta * (1 -
% delta/(beta * K)) = 0

%%

[x2, conv2, J2] = MySolve(lotkaVolterra, [0.1;0.1], lvJacobian, 10^-6, 1000);

x2

% Thus [x,y] = [0,0] is the extinction equilibrium

%% Calculate eigenvalues

eig1 = eig(J1)

% Lambda2 = 0, therefore this is a non-hyperbolic point for this
% equilibrium. A transcritical bifurcation is expected to occur at this
% equilibrium when delta / (beta * K) = 1, which the parameter values
% satisfy

%%

eig2 = eig(J2)

% This is consistent with the lambda1 = alpha and lambda2 = - delta
% determined analytically, therefore a saddle.

%% Modify parameter values

% Try delta = 0.1 and delta = 0.25 to see what happens at the non-hyperbolic
% point.

%% Delta = 0.1

delta = 0.1;

lotkaVolterra = @(x)([alpha * x(1) * (1 - x(1) / K) - beta * x(1) * x(2); beta * x(1) * x(2) - delta * x(2)]);
lvJacobian = @(x)([alpha - 2 * alpha * x(1) / K - beta * x(2), -1 * beta * x(1); beta * x(2), beta * x(1) - delta]);

[x3, conv3, J3] = MySolve(lotkaVolterra, [1.1;1.1], lvJacobian, 10^-6, 1000);

x3

eig3 = eig(J3)

% This is now only the prey-only equilibrium, and with delta=0.1, this is a
% saddle, as expected analytically, as beta * K > delta.

%%

[x4, conv4, J4] = MySolve(lotkaVolterra, [0.5;0.5], lvJacobian, 10^-6, 1000);

x4

eig4 = eig(J4)

% This is the coexistence equilibrium, and with delta=0.1, this is stable,
% as expected analytically, as delta / (beta * K) < 1.

%% Delta = 0.25

delta = 0.25;

lotkaVolterra = @(x)([alpha * x(1) * (1 - x(1) / K) - beta * x(1) * x(2); beta * x(1) * x(2) - delta * x(2)]);
lvJacobian = @(x)([alpha - 2 * alpha * x(1) / K - beta * x(2), -1 * beta * x(1); beta * x(2), beta * x(1) - delta]);

[x5, conv5, J5] = MySolve(lotkaVolterra, [1.1;1.1], lvJacobian, 10^-6, 1000);

x5

eig5 = eig(J5)

% This is now only the prey-only equilibrium, and with delta=0.25, this is a
% stable node, as expected analytically, as beta * K < delta.

%%

[x6, conv6, J6] = MySolve(lotkaVolterra, [0.5;0.5], lvJacobian, 10^-6, 1000);

x6

eig6 = eig(J6)

% This is the coexistence equilibrium, and with delta=0.25, this is a saddle,
% as expected analytically, as delta / (beta * K) > 1.

% Thus, the two equilibria exchange stability at the point delta / (beta *
% K) = 1, where there is a non-hyperbolic equilibrium, therefore the
% transcritical bifurcation has been verified.
