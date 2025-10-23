%% KYA314 - Assign 5 Q4c
% simulate the Rossler Model
% calculate Lyapunov exponents
clear;
close all;
clc;

%% define map
Rossler =@(t,x,p) [-x(2,:)-x(3,:);
                x(1,:)+p(1,:).*x(2,:);
                p(2,:)+x(3,:).*(x(1,:)-p(3,:))];

%% Parameters
% set constant parameters
b = 0.1;
c = 6;

a1 = -3; % stable equilibrium
a2 = 0.01; % stable periodic orbit
a3 = 0.2; % chaotic attractor

p = [a1,a2,a3;b,b,b;c,c,c];

%% Plot solutions
x1 = 0.1;
x2 = 0.1;
x3 = 0.1;

% Set up for initial value problem solver

x0 = [x1;x2;x3];
tspan = [0,2000];
h = 0.01;

% Solve the ODE

[X,t,xeq1] = MyIVP(@(t,x)Rossler(t,x,p),[x0,x0,x0],tspan,h);

% Plot 3D trajectory

for i =1:3
figure(1);
subplot(1,3,i); hold on;
plot3(squeeze(X(1,i,150000)),squeeze(X(2,i,150000)),squeeze(X(3,i,150000)),'k.')
plot3(squeeze(X(1,i,150000:end)),squeeze(X(2,i,150000:end)), ...
    squeeze(X(3,i,150000:end)),'k','Linewidth',2);
view(-37.5, 30);
set(gca,'FontSize',16)
xlabel('x');
ylabel('y');
zlabel('z');

% Plot trajectory in time

figure(2);
subplot(1,3,i);
plot(t(150000:end),squeeze(X(:,i,150000:end)),'Linewidth',2);
set(gca,'FontSize',16)
xlabel('t');
legend('x','y','z')
xlim([1500,2000])
end


