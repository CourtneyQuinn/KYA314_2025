%% Forced Van der Pol Oscillator
clc;
clear;
close all;

%% Define the model
% Here we model the forced Duffing oscillator
Duff=@(t,x,p)[x(2,:); % velocity
            -p(1,:).*x(2,:)-p(2,:).*x(1,:)-p(3,:).*x(1,:).^3+p(4,:).*cos(p(5,:).*t)]; % main Duffing equation

%% Set parameters
% 𝛿 = 0,    α = 1,     β = 0, K = 0.5,  ω = 0.5       P
% 𝛿 = 0,    α = 1,     β = 0, K = 0.5,  ω = 1.2       P
% 𝛿 = 0,    α = 1,     β = 0, K = 0.5,  ω = sqrt(2)   QP
% 𝛿 = 0,    α = 1,     β = 1, K = 0.5,  ω = 0.5       QP
% 𝛿 = 1,    α = 1,     β = 0, K = 0.5,  ω = 0.5       P
% 𝛿 = 1,    α = 1,     β = 1, K = 10,   ω = 1.2       Chaos

delta = 0; % damping
alpha = 1; % stiffness
beta = 1; % nonlinearity
K = 10; % forcing strength
omega = 0.5; % frequency

p = [delta;alpha;beta;K;omega];

%% Solve IVP

x1 = 0;
x2 = 2;

% Set up for initial value problem solver

x0 = [x1;x2];
tspan = [0,1000];
timescale = 2/omega*pi;
h = 0.01*timescale;

% Solve the ODE

[X,t,xeq1] = MyIVP(@(t,x)Duff(t,x,p),x0,tspan,h);

figure(1); 
plot(t,X(1,:));
set(gca,'FontSize',16)
xlabel('t');
ylabel('x');

%% Calculate frequency spectrum

fs = 1/h;
N = length(X(1,:));
xdft = fft(X(1,:));
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(X(1,:)):fs/2;

figure(4);
plot(2*pi*freq,psdx./max(psdx).*freq)
grid on
title("Spectrum")
xlabel("Frequency")
ylabel("Power spectral density")
xlim([0 5])

[~,max_loc]=find(psdx.*freq == max(psdx.*freq));
max_freq = freq(max_loc)*2*pi

%% Plot poincare map

%step = floor(max_freq/0.01);

figure(2); 
plot(X(1,1000:100:end),X(2,1000:100:end),'.','MarkerSize',8,'Linewidth',2);
set(gca,'FontSize',16)
xlabel('x');
ylabel('v');
%axis([-4 4 -4 4])

%% Plot 3D trajectory

figure(3);
plot3(X(1,100:end),X(2,100:end),cos(omega*t(100:end)),'.','MarkerSize',8,'Linewidth',2);
set(gca,'FontSize',16)
xlabel('x');
ylabel('v');
zlabel('z');

%% Calculate Lyapunov exponents
hjac = 1e-6;
Js = MyJacobian(@(x)Duff(0,x,p),X,hjac);
M = NaN(size(Js));
for j = 1:size(Js,3)
    M(:,:,j) = expm(Js(:,:,j)*h);
end

N = size(X,2)-1;

[ lambda,~,~,~] = LyapQR_new(M,x0,N,[],h);
disp(lambda)



