%% KYA314 - Assign 5 Q1c 
% Duffing oscillator
% Perturbed homoclinic orbit
clear 
clc
close all

%% Define system
f=@(t,x,p)[x(2,:);
    -p(1,:).*x(1,:)-p(2,:).*x(1,:).^3+p(3,:).*cos(t)];

alpha=-1;
beta=0.25;

%% No forcing
p=[alpha;beta;0];

% ODE solver
tspan = [0,100];
h = 1e-2;
x0 = [1e-3;1e-3];

[xout1,t,~] = MyIVP(@(t,x)f(t,x,p),x0,tspan,h);

% Plot in phase space and time on same figure
figure(1); hold on;
% plot in phase space
subplot(1,2,1)
plot(xout1(1,:),xout1(2,:))
set(gca,'FontSize',16)
xlabel('$x$','interpreter','latex')
ylabel('$\dot{x}$','interpreter','latex')
% plot in time
subplot(1,2,2); hold on;
plot(t,xout1)
set(gca,'FontSize',16)
xlabel('$t$','interpreter','latex')
legend('$x$','$\dot{x}$','interpreter','latex')
box on;

%% Weak forcing (eps = 0.01)
p=[alpha;beta;0.01];

% ODE solver
tspan = [0,100];
h = 1e-2;
x0 = [1e-3;1e-3];

[xout1,t,~] = MyIVP(@(t,x)f(t,x,p),x0,tspan,h);

% Plot in phase space and time on same figure
figure(2); hold on;
% plot in phase space
subplot(1,2,1)
plot(xout1(1,:),xout1(2,:))
set(gca,'FontSize',16)
xlabel('$x$','interpreter','latex')
ylabel('$\dot{x}$','interpreter','latex')
% plot in time
subplot(1,2,2); hold on;
plot(t,xout1)
set(gca,'FontSize',16)
xlabel('$t$','interpreter','latex')
legend('$x$','$\dot{x}$','interpreter','latex')
box on;

%% Interpretation
% The second case is experiencing a pertubation of a homocliniic orbit
% which leads to a chaotic solution.