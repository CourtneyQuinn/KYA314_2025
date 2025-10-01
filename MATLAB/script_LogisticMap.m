%% KYA314 - Logistic Map
% simulate the recursive sequence for varying values of r
% calculate Lyapunov exponents
clear;
close all;
clc;

%% r = 3.57 (chaos)
% set parameters
r = 3.57;
x0 = 0.1;
Nsteps = 1000;

% create empty solution vector 
xtraj = NaN(Nsteps+1,1);

% input initial condition
xtraj(1,:) = x0;

% iterate map
for i = 1:Nsteps
    x0 = LogisticMap(x0,r);
    xtraj(i+1,:) = x0;
end

% plot solution
figure(1); clf;
plot(xtraj(100:end-1,:),xtraj(101:end,:),'.','MarkerSize',10,'Linewidth',3)
axis([-1 1 -1 1])
xlabel("x_n")
ylabel("x_{n+1}")
title("Solution trajectory in phase space")

figure(2); clf;
plot(linspace(101,Nsteps+1,Nsteps-99),xtraj(101:end,:),'.','MarkerSize',8,'Linewidth',3)
xlim([100 1000])
xlabel("n")
ylabel("x_n")
title("Solution trajectory")

%% Calculate Lyapunov exponent
tan_func =@(x) r.*(1-2.*x);
logdf =@(x) log(abs(tan_func(x)));

Lyap = mean(logdf(xtraj));

disp("Lyap exp for r=" + num2str(r) + " is " + num2str(Lyap))

%% Plot Lyap exp for varying r

r_vals = linspace(1,4,300);

Lyaps = NaN(length(r_vals),1);

for j = 1:length(r_vals)
    r = r_vals(j);
    
    % create empty solution vector 
    xtraj = NaN(Nsteps+1,1);
    
    % input initial condition
    xtraj(1,:) = x0;
    
    % iterate map
    for i = 1:Nsteps
        x0 = LogisticMap(x0,r);
        xtraj(i+1,:) = x0;
    end
    
    tan_func =@(x) r.*(1-2.*x);
    logdf =@(x) log(abs(tan_func(x)));

    Lyap = mean(logdf(xtraj));

    Lyaps(j) = Lyap;
end

figure(3); clf;
hold on;
plot([1,4],[0,0],'k-','Linewidth',2)
plot(r_vals,Lyaps,'r.','MarkerSize',8,'Linewidth',3)
xlabel("r")
ylabel("\lambda")
title("Lyapunov exponents")
