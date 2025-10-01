%% KYA314 - Intermittent chaos
% simulate intermittent behaviour
% calculate Lyapunov exponents
clear;
close all;
clc;

%% define map
f =@(x,r) (r + x - x.^2)./(1+x.^3);

%% Check r values as approaching zero from above 0 
r = 0.1;
x01 = sqrt(-1/2+1/2*sqrt(1+4*r));
x02 = -sqrt(-1/2+1/2*sqrt(1+4*r));
Nsteps = 50;

% create empty solution vector 
xtraj = NaN(Nsteps+1,2);

% input initial condition
xtraj(1,:) = [x01;x02];

% iterate map
for i = 1:Nsteps
    x01 = f(x01,r);
    x02 = f(x02,r);
    xtraj(i+1,:) = [x01;x02];
end

% plot solution
figure(1); clf;
subplot(1,2,1)
plot(xtraj(1:end-1,:),xtraj(2:end,:),'.','MarkerSize',10,'Linewidth',3)
xlabel("x_n")
ylabel("x_{n+1}")
title("Solution trajectory in phase space")
axis([-1 1 -1 1])

subplot(1,2,2)
plot(linspace(1,Nsteps+1,Nsteps+1),xtraj,'.','MarkerSize',8,'Linewidth',3)
xlabel("n")
ylabel("x_n")
title("Solution trajectory")
ylim([-1 1])

%% Check r values below 0 
%  (check r=-1?)
r = -0.001;
x0 = 0;
Nsteps = 100;

% create empty solution vector 
xtraj = NaN(Nsteps+1,1);

% input initial condition
xtraj(1,:) = x0;

% iterate map
for i = 1:Nsteps
    x0 = f(x0,r);
    xtraj(i+1,:) = x0;
end

% plot solution
figure(2); clf;
subplot(1,2,1)
plot(xtraj(1:end-1,:),xtraj(2:end,:),'.','MarkerSize',10,'Linewidth',3)
xlabel("x_n")
ylabel("x_{n+1}")
title("Solution trajectory in phase space")

subplot(1,2,2)
plot(linspace(1,Nsteps+1,Nsteps+1),xtraj,'.','MarkerSize',8,'Linewidth',3)
xlabel("n")
ylabel("x_n")
title("Solution trajectory")

%% Calculate Lyapunov exponent
tan_func =@(x,r) (1-2.*x-3.*r.*x.^2-2.*x.^3+x.^4)./(1+x.^3).^2;
logdf =@(x,r) log(abs(tan_func(x,r)));

Lyap = mean(logdf(xtraj,r));

disp("Lyap exp for r=" + num2str(r) + " is " + num2str(Lyap))


%% Plot Lyapunov exponents for varying r
r_vals = linspace(-2,2,400);

Lyaps = NaN(length(r_vals),1);
Nsteps = 1000;

for j = 1:length(r_vals)
    r = r_vals(j);
    x0 = 0.001;

    % create empty solution vector 
    xtraj = NaN(Nsteps+1,1);
    
    % input initial condition
    xtraj(1,:) = x0;
    
    % iterate map
    for i = 1:Nsteps
        x0 = f(x0,r);
        xtraj(i+1,:) = x0;
    end
    
    tan_func =@(x) (1-2.*x-3.*r.*x.^2-2.*x.^3+x.^4)./(1+x.^3).^2;
    logdf =@(x) log(abs(tan_func(x)));

    Lyap = mean(logdf(xtraj));

    Lyaps(j) = Lyap;
end

figure(3); clf;
hold on;
plot([-2,2],[0,0],'k-','Linewidth',2)
plot(r_vals,Lyaps,'r.','MarkerSize',8,'Linewidth',3)
xlabel("r")
ylabel("\lambda")
title("Lyapunov exponents")