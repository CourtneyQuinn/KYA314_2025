%% KYA314 - Assignment 1, Q5de
% simulate the Duffing Map for given values of alpha and beta
clear;
close all;
clc;

% define map
f=@(x,p) [x(2,:);
          -p(2).*x(1,:)+p(1).*x(2,:)-x(2,:).^3];

%% Part d
% set parameters
alpha = 1.4;
beta = 0.25;
p = [alpha;beta];
Nsteps = 100;
x1 = -1;
y1 = -1;
xtemp1 = [x1;y1];
x2 = 1;
y2 = 1;
xtemp2 = [x2;y2];

% create empty solution vector 
xtraj1 = NaN(2,Nsteps+1);
xtraj2 = NaN(2,Nsteps+1);
xtraj1(:,1) = xtemp1;
xtraj2(:,1) = xtemp2;

% iterate map
for i = 1:Nsteps
    xtemp1 = f(xtemp1,p);
    xtraj1(:,i+1) = xtemp1;
    xtemp2 = f(xtemp2,p);
    xtraj2(:,i+1) = xtemp2;
end

% plot solution
figure(1); clf; hold on;
plot(xtraj1(1,50:end),xtraj1(2,50:end),'.','MarkerSize',12,'Linewidth',3)
plot(xtraj2(1,50:end),xtraj2(2,50:end),'.','MarkerSize',12,'Linewidth',3)
%axis([-2 2 -0.5 0.5])
set(gca,'FontSize',12)
xlabel('x_n')
ylabel('y_n')
title("Trajectories in Phase Plane")

figure(2); clf; hold on;
plot(linspace(0,Nsteps,Nsteps+1),xtraj1(1,:),'.','MarkerSize',8,'Linewidth',3)
plot(linspace(0,Nsteps,Nsteps+1),xtraj2(1,:),'.','MarkerSize',8,'Linewidth',3)
xlim([0 100])
set(gca,'FontSize',12)
xlabel('n')
ylabel('x_n')
title("Solution Trajectories")

%% Observation for part d)
% The orbits converge to stable points, interestingly in the opposite 
% quadrant to where they start
%% Chaos
% set parameters
alpha = 2.1;
beta = 0.25;
p = [alpha;beta];
Nsteps = 100;
x1 = -1;
y1 = -1;
xtemp1 = [x1;y1];
x2 = 1;
y2 = 1;
xtemp2 = [x2;y2];

% create empty solution vector 
xtraj1 = NaN(2,Nsteps+1);
xtraj2 = NaN(2,Nsteps+1);
xtraj1(:,1) = xtemp1;
xtraj2(:,1) = xtemp2;

% iterate map
for i = 1:Nsteps
    xtemp1 = f(xtemp1,p);
    xtraj1(:,i+1) = xtemp1;
    xtemp2 = f(xtemp2,p);
    xtraj2(:,i+1) = xtemp2;
end

% plot solution
figure(1); clf; hold on;
plot(xtraj1(1,50:end),xtraj1(2,50:end),'.','MarkerSize',12,'Linewidth',3)
plot(xtraj2(1,50:end),xtraj2(2,50:end),'.','MarkerSize',12,'Linewidth',3)
set(gca,'FontSize',12)
xlabel('x_n')
ylabel('y_n')
title("Trajectories in Phase Plane")

figure(2); clf; hold on;
plot(linspace(0,Nsteps,Nsteps+1),xtraj1(1,:),'.','MarkerSize',8,'Linewidth',3)
plot(linspace(0,Nsteps,Nsteps+1),xtraj2(1,:),'.','MarkerSize',8,'Linewidth',3)
xlim([0 100])
set(gca,'FontSize',12)
xlabel('n')
ylabel('x_n')
title("Solution Trajectories")

%% Observation for part e)
% The orbits have swapped in which quadrant they converge. They appear to
% be converging to single points, but further investigation shows they are
% actually periodic (with small amplitude). My bad - I sho=uld have given
% you a slightly larger alpha like 2.6 to better see periodic behaviour.