%% KYA314 - Linear systems
% find eigenvectors and eigenvalues of origin for linear systems
% plot some solutions
clear;
close all;
clc;

%% 2D Examples
% set some initial points for sample trajectories
% set iteration parameters
% set empty solution vectors

x01 = [-2;0];
x02 = [0;-2];
x03 = [2;0];
x04 = [0;2];

stepsize = 0.01;
tend = 10;
Ts = linspace(0,tend,tend/stepsize);

x1s = NaN(2,length(Ts));
x2s = NaN(2,length(Ts));
x3s = NaN(2,length(Ts));
x4s = NaN(2,length(Ts));

%% stable node (sink)
A = [-5,-3;-2,-2];

[V,D] = eigs(A);

v1 = V(:,1);
v2 = V(:,2);

lams =diag(D);

lam1 = lams(1)

lam2 = lams(2)
%%
i=0;
for t = 0:stepsize:tend
    i = i+1;
    x1s(:,i) = expm(A.*t)*x01;
    x2s(:,i) = expm(A.*t)*x02;
    x3s(:,i) = expm(A.*t)*x03;
    x4s(:,i) = expm(A.*t)*x04;
end

figure(1); clf; hold on;
% eigenvectors
plot(6*exp(lam1.*Ts).*v1(1),6*exp(lam1.*Ts).*v1(2),'b','LineWidth', 1.5)
plot(-6*exp(lam1.*Ts).*v1(1),-6*exp(lam1.*Ts).*v1(2),'b','LineWidth', 1.5)

plot(6*exp(lam2.*Ts).*v2(1),6*exp(lam2.*Ts).*v2(2),'b','LineWidth', 1.5)
plot(-6*exp(lam2.*Ts).*v2(1),-6*exp(lam2.*Ts).*v2(2),'b','LineWidth', 1.5)

% example trajectories
plot(x1s(1,:),x1s(2,:),'k','LineWidth', 1.5)
plot(x2s(1,:),x2s(2,:),'k','LineWidth', 1.5)
plot(x3s(1,:),x3s(2,:),'k','LineWidth', 1.5)
plot(x4s(1,:),x4s(2,:),'k','LineWidth', 1.5)
axis([-5 5 -5 5])
xlabel('x_1')
ylabel('x_2')

%% unstable node (source)
clc;
A = [2,1;1,3];

[V,D] = eigs(A);

v1 = V(:,1);
v2 = V(:,2);

lams =diag(D);

lam1 = lams(1)

lam2 = lams(2)
%%
i=0;
for t = 0:stepsize:tend
    i = i+1;
    x1s(:,i) = expm(A.*t)*x01;
    x2s(:,i) = expm(A.*t)*x02;
    x3s(:,i) = expm(A.*t)*x03;
    x4s(:,i) = expm(A.*t)*x04;
end

figure(1); clf; hold on;
% eigenvectors
plot(1e-4*exp(lam1.*Ts).*v1(1),1e-4*exp(lam1.*Ts).*v1(2),'r','LineWidth', 1.5)
plot(-1e-4*exp(lam1.*Ts).*v1(1),-1e-4*exp(lam1.*Ts).*v1(2),'r','LineWidth', 1.5)

plot(1e-4*exp(lam2.*Ts).*v2(1),1e-4*exp(lam2.*Ts).*v2(2),'r','LineWidth', 1.5)
plot(-1e-4*exp(lam2.*Ts).*v2(1),-1e-4*exp(lam2.*Ts).*v2(2),'r','LineWidth', 1.5)

% example trajectories
plot(x1s(1,:),x1s(2,:),'k','LineWidth', 1.5)
plot(x2s(1,:),x2s(2,:),'k','LineWidth', 1.5)
plot(x3s(1,:),x3s(2,:),'k','LineWidth', 1.5)
plot(x4s(1,:),x4s(2,:),'k','LineWidth', 1.5)
axis([-5 5 -5 5])
xlabel('x_1')
ylabel('x_2')

%% saddle
clc;
A = [-1,3;2,-1];

[V,D] = eigs(A);

v1 = V(:,1);
v2 = V(:,2);

lams =diag(D);

lam1 = lams(1)

lam2 = lams(2)

%%
i=0;
for t = 0:stepsize:tend
    i = i+1;
    x1s(:,i) = expm(A.*t)*x01;
    x2s(:,i) = expm(A.*t)*x02;
    x3s(:,i) = expm(A.*t)*x03;
    x4s(:,i) = expm(A.*t)*x04;
end

figure(1); clf; hold on;
% eigenvectors
plot(10*exp(lam1.*Ts).*v1(1),10*exp(lam1.*Ts).*v1(2),'b','LineWidth', 1.5)
plot(-10*exp(lam1.*Ts).*v1(1),-10*exp(lam1.*Ts).*v1(2),'b','LineWidth', 1.5)

plot(1e-4*exp(lam2.*Ts).*v2(1),1e-4*exp(lam2.*Ts).*v2(2),'r','LineWidth', 1.5)
plot(-1e-4*exp(lam2.*Ts).*v2(1),-1e-4*exp(lam2.*Ts).*v2(2),'r','LineWidth', 1.5)

% example trajectories
plot(x1s(1,:),x1s(2,:),'k','LineWidth', 1.5)
plot(x2s(1,:),x2s(2,:),'k','LineWidth', 1.5)
plot(x3s(1,:),x3s(2,:),'k','LineWidth', 1.5)
plot(x4s(1,:),x4s(2,:),'k','LineWidth', 1.5)
axis([-5 5 -5 5])
xlabel('x_1')
ylabel('x_2')

%% stable focus
clc;
A = [-2,2;-2,-2];

[V,D] = eigs(A);

v1 = V(:,1);
v2 = V(:,2);

lams =diag(D);

lam1 = lams(1)

lam2 = lams(2)

%%
i=0;
for t = 0:stepsize:tend
    i = i+1;
    x1s(:,i) = expm(A.*t)*x01;
    x2s(:,i) = expm(A.*t)*x02;
    x3s(:,i) = expm(A.*t)*x03;
    x4s(:,i) = expm(A.*t)*x04;
end

figure(1); clf; hold on;
% eigenvectors
plot(5*exp(lam1.*Ts).*v1(1) + 5*exp(lam2.*Ts).*v2(1),...
    5*exp(lam1.*Ts).*v1(2) + 5*exp(lam2.*Ts).*v2(2),'b','LineWidth', 1.5)
plot(-5*exp(lam1.*Ts).*v1(1) - 5*exp(lam2.*Ts).*v2(1),...
    -5*exp(lam1.*Ts).*v1(2) - 5*exp(lam2.*Ts).*v2(2),'b','LineWidth', 1.5)

% example trajectories
plot(x1s(1,:),x1s(2,:),'k','LineWidth', 1.5)
plot(x2s(1,:),x2s(2,:),'k','LineWidth', 1.5)
plot(x3s(1,:),x3s(2,:),'k','LineWidth', 1.5)
plot(x4s(1,:),x4s(2,:),'k','LineWidth', 1.5)
axis([-5 5 -5 5])
xlabel('x_1')
ylabel('x_2')

%% unstable focus
clc;
A = [2,2;-2,2];

[V,D] = eigs(A);

v1 = V(:,1);
v2 = V(:,2);

lams =diag(D);

lam1 = lams(1)
lam2 = lams(2)
%%
i=0;
for t = 0:stepsize:tend
    i = i+1;
    x1s(:,i) = expm(A.*t)*x01;
    x2s(:,i) = expm(A.*t)*x02;
    x3s(:,i) = expm(A.*t)*x03;
    x4s(:,i) = expm(A.*t)*x04;
end

figure(1); clf; hold on;
% eigenvectors
plot(1e-5*exp(lam1.*Ts).*v1(1) + 1e-5*exp(lam2.*Ts).*v2(1),...
    1e-5*exp(lam1.*Ts).*v1(2) + 1e-5*exp(lam2.*Ts).*v2(2),'r','LineWidth', 1.5)
plot(-1e-5*exp(lam1.*Ts).*v1(1) - 1e-5*exp(lam2.*Ts).*v2(1),...
    -1e-5*exp(lam1.*Ts).*v1(2) - 1e-5*exp(lam2.*Ts).*v2(2),'r','LineWidth', 1.5)

% example trajectories
plot(x1s(1,:),x1s(2,:),'k','LineWidth', 1.5)
plot(x2s(1,:),x2s(2,:),'k','LineWidth', 1.5)
plot(x3s(1,:),x3s(2,:),'k','LineWidth', 1.5)
plot(x4s(1,:),x4s(2,:),'k','LineWidth', 1.5)
axis([-5 5 -5 5])
xlabel('x_1')
ylabel('x_2')

%% centre
clc;
A = [-1,3;-3,1];

[V,D] = eigs(A);

v1 = V(:,1);
v2 = V(:,2);

lams =diag(D);

lam1 = lams(1)
lam2 = lams(2)

%%
i=0;
for t = 0:stepsize:tend
    i = i+1;
    x1s(:,i) = expm(A.*t)*x01;
    x2s(:,i) = expm(A.*t)*x02;
    x3s(:,i) = expm(A.*t)*x03;
    x4s(:,i) = expm(A.*t)*x04;
end

figure(1); clf; hold on;
% eigenvectors
plot(1e-1*exp(lam1.*Ts).*v1(1) + 1e-1*exp(lam2.*Ts).*v2(1),...
    1e-1*exp(lam1.*Ts).*v1(2) + 1e-1*exp(lam2.*Ts).*v2(2),'g','LineWidth', 1.5)
plot(exp(lam1.*Ts).*v1(1) + exp(lam2.*Ts).*v2(1),...
    exp(lam1.*Ts).*v1(2) + exp(lam2.*Ts).*v2(2),'g','LineWidth', 1.5)

% example trajectories
plot(x1s(1,:),x1s(2,:),'k','LineWidth', 1.5)
plot(x2s(1,:),x2s(2,:),'k','LineWidth', 1.5)
plot(x3s(1,:),x3s(2,:),'k','LineWidth', 1.5)
plot(x4s(1,:),x4s(2,:),'k','LineWidth', 1.5)
axis([-5 5 -5 5])
xlabel('x_1')
ylabel('x_2')

%% 3D Examples
% create grid in x for plotting
% set some initial points
% set empty solution vectors

x01 = [2;2;0];
x02 = [0;2;1];
x03 = [2;0;2];

x1s = NaN(3,length(Ts));
x2s = NaN(3,length(Ts));
x3s = NaN(3,length(Ts));

%% 3D Example 1
clc;
A = [-1,-2,0;1,-2,-2;0,1,-2];

[V,D] = eigs(A);

v1 = V(:,1);
v2 = V(:,2);
v3 = V(:,3);

lams =diag(D);

lam1 = lams(1)
lam2 = lams(2)
lam3 = lams(3)

%%
i=0;
for t = 0:stepsize:tend
    i = i+1;
    x1s(:,i) = expm(A.*t)*x01;
    x2s(:,i) = expm(A.*t)*x02;
    x3s(:,i) = expm(A.*t)*x03;
end

figure(2); clf; hold on;
% eigenvectors
plot3(exp(lam1.*Ts).*v1(1) + exp(lam2.*Ts).*v2(1),...
    exp(lam1.*Ts).*v1(2) + exp(lam2.*Ts).*v2(2),...
    exp(lam1.*Ts).*v1(3) + exp(lam2.*Ts).*v2(3),'b','LineWidth', 1.5)
plot3(-exp(lam1.*Ts).*v1(1) - exp(lam2.*Ts).*v2(1),...
    -exp(lam1.*Ts).*v1(2) - exp(lam2.*Ts).*v2(2),...
    -exp(lam1.*Ts).*v1(3) - exp(lam2.*Ts).*v2(3),'b','LineWidth', 1.5)

plot3(exp(lam3.*Ts).*v3(1),...
    exp(lam3.*Ts).*v3(2),...
    exp(lam3.*Ts).*v3(3),'b','LineWidth', 1.5)
plot3(-exp(lam3.*Ts).*v3(1),...
    -exp(lam3.*Ts).*v3(2),...
    -exp(lam3.*Ts).*v3(3),'b','LineWidth', 1.5)


% example trajectories
plot3(x1s(1,:),x1s(2,:),x1s(3,:),'k','LineWidth', 1.5)
plot3(x2s(1,:),x2s(2,:),x2s(3,:),'k','LineWidth', 1.5)
plot3(x3s(1,:),x3s(2,:),x3s(3,:),'k','LineWidth', 1.5)
axis([-1 1 -1 1 -1 1])
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')

%% 3D Example 2
clc;
A = [0,-2,0;2,0,-2;0,2,0];

[V,D] = eigs(A);

v1 = V(:,1);
v2 = V(:,2);
v3 = V(:,3);

lams =diag(D);

lam1 = lams(1)
lam2 = lams(2)
lam3 = lams(3)

%%
i=0;
for t = 0:stepsize:tend
    i = i+1;
    x1s(:,i) = expm(A.*t)*x01;
    x2s(:,i) = expm(A.*t)*x02;
    x3s(:,i) = expm(A.*t)*x03;
end

figure(2); clf; hold on;
% eigenvectors
plot3(exp(lam1.*Ts).*v1(1) + exp(lam2.*Ts).*v2(1),...
    exp(lam1.*Ts).*v1(2) + exp(lam2.*Ts).*v2(2),...
    exp(lam1.*Ts).*v1(3) + exp(lam2.*Ts).*v2(3),'b','LineWidth', 1.5)

plot3(exp(lam3.*Ts).*v3(1),...
    exp(lam3.*Ts).*v3(2),...
    exp(lam3.*Ts).*v3(3),'b','LineWidth', 1.5)

% example trajectories
plot3(x1s(1,:),x1s(2,:),x1s(3,:),'k','LineWidth', 1.5)
plot3(x2s(1,:),x2s(2,:),x2s(3,:),'k','LineWidth', 1.5)
plot3(x3s(1,:),x3s(2,:),x3s(3,:),'k','LineWidth', 1.5)
axis([-4 4 -4 4 -4 4])
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')