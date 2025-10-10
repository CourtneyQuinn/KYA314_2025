%% KYA314 Assignment 4, Q4

clear;
clc;
format long;

%% Part A

disp("Part A")

% Define parameters

K = 6;
delta = 1.3;

% Define the predator-prey system

predatorPrey = @(t,x)([x(1) * (1 - x(1) / K) - x(1) * x(2) / (1 + 0.1 * x(1) * x(1)); x(1) * x(2) / (1 + 0.1 * x(1) * x(1)) - delta * x(2)]);

% Run MySolve to find the equilibria

[xt, t, xend] = MyIVP(predatorPrey, [1;1], [0,150], 1e-3);

figure(1);
hold on;
plot(t, xt(1,:));
plot(t, xt(2,:));
xlabel('t')
legend('x', 'y')
title("Solution trajectory in time")
hold off;

figure(2);
hold on;
plot(xt(1,:), xt(2,:));
xlabel('x')
ylabel('y')
title("Solution trajectory in phase space")
hold off;

disp("This appears to be an unstable focus converging to a periodic orbit.")

%% Part B

disp("Part B")

% Find the coexistence equilibrium

predatorPreyNoT = @(x)([x(1) * (1 - x(1) / K) - x(1) * x(2) / (1 + 0.1 * x(1) * x(1)); x(1) * x(2) / (1 + 0.1 * x(1) * x(1)) - delta * x(2)]);

df=@(x)MyJacobian(predatorPreyNoT,x,1e-3);

[x1,~] = MySolve(predatorPreyNoT,xend,df,1e-6,1000);

disp("The equilibrium is x = " + x1(1) + ", y = " + x1(2));

J1 = df(x1);

eig1 = eig(J1);

disp("The eigenvalues of the equilibrium are " + eig1(1) + ", " + eig1(2));

figure(3);
hold on;
plot(xt(1,:), xt(2,:));
plot(x1(1), x1(2), "*");
xlabel('x')
ylabel('y')
title("Solution trajectory in phase space")
legend("Solution trajectory", "Equilibrium")
hold off;

disp("Thus the coexistence equilibrium is an unstable focus.")

%% Part C

disp("Part C")

K = 5.5;
delta = 1.3;

predatorPrey = @(t,x)([x(1) * (1 - x(1) / K) - x(1) * x(2) / (1 + 0.1 * x(1) * x(1)); x(1) * x(2) / (1 + 0.1 * x(1) * x(1)) - delta * x(2)]);

[xt, t, xend] = MyIVP(predatorPrey, x1, [0,150], 1e-3);

figure(4);
hold on;
plot(t, xt(1,:));
plot(t, xt(2,:));
xlabel('t')
legend('x', 'y')
title("Solution trajectory in time")
hold off;

figure(5);
hold on;
plot(xt(1,:), xt(2,:));
xlabel('x')
ylabel('y')
title("Solution trajectory in phase space")
hold off;

predatorPreyNoT = @(x)([x(1) * (1 - x(1) / K) - x(1) * x(2) / (1 + 0.1 * x(1) * x(1)); x(1) * x(2) / (1 + 0.1 * x(1) * x(1)) - delta * x(2)]);

df=@(x)MyJacobian(predatorPreyNoT,x,1e-3);

[x1,~] = MySolve(predatorPreyNoT,xend,df,1e-6,1000);

disp("The equilibrium is x = " + x1(1) + ", y = " + x1(2));

J1 = df(x1);

eig1 = eig(J1);

disp("The eigenvalues of the equilibrium are " + eig1(1) + ", " + eig1(2));

figure(6);
hold on;
plot(xt(1,:), xt(2,:));
plot(x1(1), x1(2), "*");
xlabel('x')
ylabel('y')
title("Solution trajectory in phase space")
legend("Solution trajectory", "Equilibrium")
hold off;

disp("This time the coexistence equilibrium appears to be a stable focus with very slow convergence.")
disp("The original periodic orbit no longer exists.")
disp("Therefore we suspect a supercritical Hopf bifurcation.")

%% Part D

disp("Part D")

K = 5.5;
delta = 1 / (2 * sqrt(0.1));

predatorPrey = @(t,x)([x(1) * (1 - x(1) / K) - x(1) * x(2) / (1 + 0.1 * x(1) * x(1)); x(1) * x(2) / (1 + 0.1 * x(1) * x(1)) - delta * x(2)]);

[xt, t, xend] = MyIVP(predatorPrey, x1, [0,150], 1e-3);

figure(7);
hold on;
plot(t, xt(1,:));
plot(t, xt(2,:));
xlabel('t')
legend('x', 'y')
title("Solution trajectory in time")
hold off;

figure(8);
hold on;
plot(xt(1,:), xt(2,:));
xlabel('x')
ylabel('y')
title("Solution trajectory in phase space")
hold off;

predatorPreyNoT = @(x)([x(1) * (1 - x(1) / K) - x(1) * x(2) / (1 + 0.1 * x(1) * x(1)); x(1) * x(2) / (1 + 0.1 * x(1) * x(1)) - delta * x(2)]);

df=@(x)MyJacobian(predatorPreyNoT,x,1e-3);

[x1,~] = MySolve(predatorPreyNoT,x1,df,1e-6,1000);

disp("The equilibrium is x = " + x1(1) + ", y = " + x1(2));

J1 = df(x1);

eig1 = eig(J1);

disp("The eigenvalues of the equilibrium are " + eig1(1) + ", " + eig1(2));

figure(9);
hold on;
plot(xt(1,:), xt(2,:));
plot(x1(1), x1(2), "*");
xlabel('x')
ylabel('y')
title("Solution trajectory in phase space")
legend("Solution trajectory", "Equilibrium")
hold off;

disp("The coexistence equilibrium is experiencing a fold (saddle-node) bifurcation.")

%% Part E

disp("Part E")

K = 2 / sqrt(0.1);
delta = 1 / (2 * sqrt(0.1));

predatorPrey = @(t,x)([x(1) * (1 - x(1) / K) - x(1) * x(2) / (1 + 0.1 * x(1) * x(1)); x(1) * x(2) / (1 + 0.1 * x(1) * x(1)) - delta * x(2)]);

[xt, t, xend] = MyIVP(predatorPrey, x1, [0,150], 1e-3);

figure(10);
hold on;
plot(t, xt(1,:));
plot(t, xt(2,:));
xlabel('t')
legend('x', 'y')
title("Solution trajectory in time")
hold off;

figure(11);
hold on;
plot(xt(1,:), xt(2,:));
xlabel('x')
ylabel('y')
title("Solution trajectory in phase space")
hold off;

predatorPreyNoT = @(x)([x(1) * (1 - x(1) / K) - x(1) * x(2) / (1 + 0.1 * x(1) * x(1)); x(1) * x(2) / (1 + 0.1 * x(1) * x(1)) - delta * x(2)]);

df=@(x)MyJacobian(predatorPreyNoT,x,1e-3);

[x1,~] = MySolve(predatorPreyNoT,x1,df,1e-6,1000);

disp("The equilibrium is x = " + x1(1) + ", y = " + x1(2));

J1 = df(x1);

eig1 = eig(J1);

disp("The eigenvalues of the equilibrium are " + eig1(1) + ", " + eig1(2));

figure(12);
hold on;
plot(xt(1,:), xt(2,:));
plot(x1(1), x1(2), "*");
xlabel('x')
ylabel('y')
title("Solution trajectory in phase space")
legend("Solution trajectory", "Equilibrium")
hold off;

disp("The coexistence equilibrium is experiencing a Bogdanov-Takens (double-zero eigenvalue) bifurcation.")