%% VanderPol oscillator
% visualise relaxation oscillations
clear;
close all;
clc;

% define ODE system
f=@(t,x,p) [x(2,:);
          p(1).*(1-x(1,:).^2).*x(2,:)-x(1,:)];

%% alpha < 0
alpha = -1;

x0 = [0.1;0.1];

% ODE solver
tspan = [0,100];
h = 1e-3;

[xout1,t,~] = MyIVP(@(t,x)f(t,x,alpha),x0,tspan,h);

% Plot in phase space and time on same figure
figure(1); clf;
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

%% alpha = 0
alpha = 0;

x0 = [0.1;0.1];

% ODE solver
tspan = [0,100];

[xout2,t,~] = MyIVP(@(t,x)f(t,x,alpha),x0,tspan,h);

% Plot in phase space
figure(1); clf;
% plot in phase space
subplot(1,2,1)
plot(xout2(1,:),xout2(2,:))
set(gca,'FontSize',16)
xlabel('$x$','interpreter','latex')
ylabel('$\dot{x}$','interpreter','latex')
% plot in time
subplot(1,2,2); hold on;
plot(t,xout2)
set(gca,'FontSize',16)
xlabel('$t$','interpreter','latex')
legend('$x$','$\dot{x}$','interpreter','latex')
box on;

%% alpha > 0
alpha = 1;

x0 = [0.1;0.1];

% ODE solver
tspan = [0,100];

[xout3,t,~] = MyIVP(@(t,x)f(t,x,alpha),x0,tspan,h);

% Plot in phase space
figure(1); clf;
% plot in phase space
subplot(1,2,1)
plot(xout3(1,:),xout3(2,:))
set(gca,'FontSize',16)
xlabel('$x$','interpreter','latex')
ylabel('$\dot{x}$','interpreter','latex')
% plot in time
subplot(1,2,2); hold on;
plot(t,xout3)
set(gca,'FontSize',16)
xlabel('$t$','interpreter','latex')
legend('$x$','$\dot{x}$','interpreter','latex')
box on;

%% alpha >> 0
alpha = 10;

x0 = [0.1;0.1];

% ODE solver
tspan = [0,100];

[xout4,t,~] = MyIVP(@(t,x)f(t,x,alpha),x0,tspan,h);

% Plot in phase space
figure(1); clf;
% plot in phase space
subplot(1,2,1)
plot(xout4(1,:),xout4(2,:))
set(gca,'FontSize',16)
xlabel('$x$','interpreter','latex')
ylabel('$\dot{x}$','interpreter','latex')
% plot in time
subplot(1,2,2); hold on;
plot(t,xout4)
set(gca,'FontSize',16)
xlabel('$t$','interpreter','latex')
legend('$x$','$\dot{x}$','interpreter','latex')
box on;
