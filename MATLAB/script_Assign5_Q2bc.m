%% KYA314 - Assign 5 Q2bc
% Reduced FitzHugh Nagumo oscillator
% visualise relaxation oscillations
clear;
close all;
clc;

%% define ODE system
f=@(t,x,p) [(1-x(1,:).^2./3).*x(1,:)-x(2,:);
          p(2,:).*(x(1,:)+p(1,:))];

%% b-c)
alpha = 0.9;

ps=[alpha,alpha,alpha;1,0.1,0.01];

x0 = [-0.1;0.1];

% Calculate slow manifold
v_slow = linspace(-2,2,100);
w = v_slow-v_slow.^3./3;

tspan = [0,1000];
h = 1e-3;

[xout,t,~] = MyIVP(@(t,x)f(t,x,ps),[x0,x0,x0],tspan,h);

for i =1:3
% Plot in phase space and time on same figure
figure(i);
% plot in phase space
subplot(1,2,1); hold on;
plot(v_slow,w,'k')
plot(squeeze(xout(1,i,:)),squeeze(xout(2,i,:)))
set(gca,'FontSize',16)
xlabel('$v$','interpreter','latex')
ylabel('$w$','interpreter','latex')
% plot in time
subplot(1,2,2); hold on;
plot(t,squeeze(xout(:,i,:)))
set(gca,'FontSize',16)
xlabel('$t$','interpreter','latex')
legend('$v$','$w$','interpreter','latex')
box on;
end

%% Interpretation
% As epsilon is decreased, the system spends more time along the slow
% manifold and quickly transitions through phse space when the slow
% manifold has a turning point. The behaviour is a relaxation oscillation
% occurring on two different timescales.