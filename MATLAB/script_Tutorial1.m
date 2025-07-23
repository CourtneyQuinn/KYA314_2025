%% KYA314 - Logistic Map
% simulate the discrete dynamical system for varying values of lam
clear;
close all;
clc;

%% set parameters
lam = 1.8;
x0 = 0.1;
Nsteps = 1000;

% create empty solution vector 
xtraj = NaN(Nsteps+1,1);

% input initial condition
xtraj(1,:) = x0;

% iterate map
for i = 1:Nsteps
    x0 = LogisticMap(x0,lam);
    xtraj(i+1,:) = x0;
end

%% plot solution

% orbit in phase space
figure(1); clf;
plot(xtraj(100:end-1,:),xtraj(101:end,:),'.','MarkerSize',10,'Linewidth',3)
axis([-1 1 -1 1])
xlabel("x_n")
ylabel("x_{n+1}")
title("Solution trajectory in phase space")

% orbit in time
figure(2); clf;
plot(linspace(101,Nsteps+1,Nsteps-99),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
xlim([100 1000])
xlabel("n")
ylabel("x_n")
title("Solution trajectory")
