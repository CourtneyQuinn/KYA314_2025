%% Script to test LyapQR
clear;
clc;

%% Henon Map

Henon =@(x,p)[1-p(1,:).*x(1,:)^2+x(2,:);
             p(2,:).*x(1,:)];

%% Iterate map for expected period-2 solution
x0 = [1;1];
p = [0.8;0.4];
N = 1000;

tn = linspace(1,N,N);
xn = NaN(size(x0,1),size(tn,2));
xn(:,1) = x0;

for ti = tn(1:N)
    xi = Henon(xn(:,ti),p);
    xn(:,ti+1) = xi;
end

plot(xn(1,end-10:end),xn(2,end-10:end),'.','Markersize',12);

%% Calculate Lyapunov exponents
[lambda,Rdiag,x,~] = LyapQR_new(@(x)Henon(x,p),x0,N);

disp(exp(2*lambda))

%% Check time-2 map 
M2 =@(x)Henon(Henon(x,p),p);

J = MyJacobian(M2,xn(:,end),1e-6);

[V,D] = eigs(J);

disp(abs(diag(D)))

%% Change to alpha = 1.4 & beta = 0.3
x0 = [1;1];
p = [1.4;0.3];
N = 1000;

tn = linspace(1,N,N);
xn = NaN(size(x0,1),size(tn,2));
xn(:,1) = x0;

for ti = tn(1:N)
    xi = Henon(xn(:,ti),p);
    xn(:,ti+1) = xi;
end

plot(xn(1,100:end),xn(2,100:end),'.');

%% Calculate Lyapunov exponents
[lambda,Rdiag,x,~] = LyapQR_new(@(x)Henon(x,p),x0,N);

disp(lambda)