%% KYA314 - Assignment 4, Q2
% plot Hamiltonian curve for unforced simple pendulum
clear;
close all;
clc;

%% Hamiltonian given by
% $ H(\theta,v) = -\frac{v^2}{2}+\alpha^2\cos(\theta).$
% 
% At saddle, $H=\alpha^2$ 
% 
% We can rearrange for $v$:
% 
% $v(\theta) = \pm\sqrt{2\alpha^2\cos(\theta)+2\alpha^2}.$
%
alpha1 = 0.5;
alpha2 = 1;
alpha3 = 2;

% Define v as function of theta
v =@(theta,alpha) sqrt(2*alpha^2*cos(theta)+2*alpha^2);

thetas = linspace(-pi,pi,628);
v1 = v(thetas,alpha1);
v2 = v(thetas,alpha2);
v3 = v(thetas,alpha3);

%% Plot
figure(1); clf;
hold on;
plot(thetas,v1)
plot(thetas,v2)
plot(thetas,v3)
plot(thetas,-v1,'Color',"#0072BD")
plot(thetas,-v2,'Color',"#D95319")
plot(thetas,-v3,'Color',"#EDB120")
xlabel("\theta")
ylabel("v")
title("Hamiltonians in Phase Plane")
xlim([-pi, pi])
ylim([-5,5])
legend('\alpha=0.5','\alpha=1','\alpha=2')
