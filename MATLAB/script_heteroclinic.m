%% Heteroclinic connection example
% This script calculates a heteroclinic orbit and its destruction using a
% forward and backward integrator
clear;
close all;
clc;

%% Define the model
% from lectures

f=@(t,x,delta)[x(2,:);
    x(1,:).^3-x(1,:)+delta.*x(2,:)];

%% delta=0 case
delta = 0;
hjac = 1e-6;

% find eigenvectors of saddles
x1 = [-1;0];
J1 = MyJacobian(@(x)f(0,x,delta),x1,hjac);
[V1,D1] = eigs(J1);

x2 = [1;0];
J2 = MyJacobian(@(x)f(0,x,delta),x2,hjac);
[V2,D2] = eigs(J2);

%% forward integrate from saddle (-1,0) along unstable eigenvector
tspan = [0,20];
h = 1e-2;
[xt1,t1,~] = MyIVP(@(t,x)f(t,x,delta),x1+V1(:,1)*1e-12,tspan,h);

% forward integrate from saddle (1,0) along unstable eigenvector
tspan = [0,20];
h = 1e-2;
[xt2,t2,~] = MyIVP(@(t,x)f(t,x,delta),x2-V2(:,1)*1e-12,tspan,h);

% plot behaviour in phase space
figure(1); clf; hold on;
plot(xt1(1,:),xt1(2,:))
plot(xt2(1,:),xt2(2,:))
axis([-1 1 -1 1])

% plot behaviour in time
figure(2); clf; hold on;
plot(t1,xt1(1,:))


%% delta>0 case
delta = 0.1;
hjac = 1e-6;

% find eigenvectors of saddles
x1 = [-1;0];
J1 = MyJacobian(@(x)f(0,x,delta),x1,hjac);
[V1,D1] = eigs(J1);

x2 = [1;0];
J2 = MyJacobian(@(x)f(0,x,delta),x2,hjac);
[V2,D2] = eigs(J2);

%% forward integrate from saddle (-1,0) along unstable eigenvector
tspan = [0,20];
h = 1e-2;
[xt1,t1,~] = MyIVP(@(t,x)f(t,x,delta),x1-V1(:,1)*1e-8,tspan,h);

% forward integrate from saddle (1,0) along unstable eigenvector
tspan = [0,20];
h = 1e-2;
[xt2,t2,~] = MyIVP(@(t,x)f(t,x,delta),x2+V2(:,1)*1e-8,tspan,h);

figure(1); clf; hold on;
plot(xt1(1,:),xt1(2,:))
plot(xt2(1,:),xt2(2,:))
axis([-1 1 -1 1])

%% delta<0 case
delta = -0.1;
hjac = 1e-6;

% find eigenvectors of saddles
x1 = [-1;0];
J1 = MyJacobian(@(x)f(0,x,delta),x1,hjac);
[V1,D1] = eigs(J1);

x2 = [1;0];
J2 = MyJacobian(@(x)f(0,x,delta),x2,hjac);
[V2,D2] = eigs(J2);

%% forward integrate from saddle (-1,0) along unstable eigenvector
tspan = [0,100];
h = 1e-2;
[xt1,t1,~] = MyIVP(@(t,x)f(t,x,delta),x1+V1(:,2)*1e-8,tspan,h);

% forward integrate from saddle (1,0) along unstable eigenvector
tspan = [0,100];
h = 1e-2;
[xt2,t2,~] = MyIVP(@(t,x)f(t,x,delta),x2-V2(:,2)*1e-8,tspan,h);

figure(1); clf; hold on;
plot(xt1(1,:),xt1(2,:))
plot(xt2(1,:),xt2(2,:))
axis([-1 1 -1 1])