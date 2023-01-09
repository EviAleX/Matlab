% Praca domowa, zadanie 1. Meeting III
clc; % Clear command window 
clear; % Remove items from workspace, freeing up system memory

n = 100; % n. of observations
N = 10000; % n. of simulations Monte Carlo

R2_1 = zeros(N,1); % empty vector for R2 for Scenario I
R2_2 = zeros(N,1); % -/- Scenario II
t_stat_1 = zeros(N,1); % empty vector for t-statistics for Scenario I
t_stat_2 = zeros(N,1); % -/- Scenario II

%% Monte Carlo simulations
% Random walk process
for j = 1:N
    x1 = zeros(n,1); % empty vector for x1
    y1 = zeros(n,1); % empty vector for y1
    eps1 = randn(n,1); % rescale variance of eps 
    eps2 = randn(n,1); % -/-
% Generating random walk series
    i=2;
    while i<=n
        x1(i) = x1(i-1) + eps1(i);
        i = i + 1; 
    end
% Generating second random walk series
    i=2;
    while i<=n
        y1(i) = y1(i-1) + eps2(i);
        i = i + 1; 
    end
    R2_1(j,1) = calculateR2_mazur(x1, y1); % R^2 for random walk process
    t_stat_1(j,1) = calculate_t_stat_mazur(x1, y1); % t-statistics for random walk process

% White noise process 
    x2 = zeros(n,1); % empty vector for x2 
    y2 = zeros(n,1); % empty vector for y2 
    eps3 = randn(n,1); % rescale variance of eps
    eps4 = randn(n,1); % -/-
% Generating white noise series 
    i=2;
    while i<=n
        x2(i) = eps3(i);
        i = i + 1; 
    end
% Generating second white noise series 
    i=2;
    while i<=n
        y2(i) = eps4(i);
        i = i + 1; 
    end
    R2_2(j,1)= calculateR2_mazur(x2, y2); % R^2 for white noise process 
    t_stat_2(j,1)= calculate_t_stat_mazur(x2, y2); % t-statistics for white noise process 
end

%% Making plots 
% Plot #1 
subplot(2,2,1) 
histogram(R2_1, 50, Normalization="pdf"); % histogram for R2 for Scenario I, 50 - nbins 
hold on

x_1 = 0:0.001:1;
a = 1/2;
b = 49;
f_1 = x_1.^(a-1).*(1-x_1).^(b-1)/((gamma(a)*gamma(b))/gamma(a+b)); % creating beta distribution density plot, where a = (k-1)/2, b = (n-k)/2
plot(x_1, f_1, 'LineWidth', 2)
axis([0 1 0 45])

hold off
xlabel('R squared', 'FontSize', 10);
ylabel('Density', 'FontSize', 10);
legend({'histogram density estimate', 'density of \beta distribution'}, 'location', 'northeast', 'FontSize', 8);

% Plot #2 
subplot(2,2,2)
histogram(t_stat_1, Normalization="pdf") % creating histogram for t-ratios from regressing one independent random-walk 
hold on

ni1 = 100;
x_2 = -15:0.01:15; % x which will estimate the density of the t-student distribution
f_2 = (gamma((ni1+1)/2))/(sqrt(ni1*pi)*gamma(ni1/2))*(1+(x_2.^2)/ni1).^(-(ni1+1)/2); % estimate the density of the t-student distribution which has a given degree of freedom
plot(x_2, f_2,'LineWidth',2)
axis([-15 15 0 0.4])

hold off
xlabel('t-statistics', 'FontSize', 10);
ylabel('Density', 'FontSize', 10)
legend({'histogram density estimate', 'density of t-distribution'}, 'location', 'northeast', 'FontSize',7);

% Plot #3
subplot(2,2,3) 
histogram(R2_2, 30, Normalization="pdf") % histogram for white noise for Scenario I, 30 - nbins 
hold on

plot(x_1, f_1, 'LineWidth', 2)
axis([0 1 0 60])
hold off

xlabel('R squared', 'FontSize', 10);
ylabel('Density', 'FontSize', 10)
legend({'histogram density estimate', 'density of \beta distribution'}, 'location', 'northeast', 'FontSize',7);

% Plot #4
subplot(2,2,4)
histogram(t_stat_2, Normalization="pdf") % creating a histogram for t-ratios from regressing one independent white noise variable
hold on 

plot(x_2, f_2,'LineWidth',2)
axis([-15 15 0 0.4])
hold off

xlabel('t-statistics', 'FontSize', 10);
ylabel('Density', 'FontSize', 10)
legend({'histogram density estimate', 'density of t-distribution'}, 'location', 'northeast', 'FontSize',7);

%% Answers for questions 
R2_1 = sort(R2_1);
R2_2 = sort(R2_2);
k = N;
while k>0
    if R2_1(k)>0.5
        k = k-1;
    else c = k;
        k = 0;
    end
end
p_1 = 1-c/N;
k = N;
while k>0
    if R2_2(k)>0.5
        k = k-1;
    else c = k;
        k = 0;
    end
end
p_2 = 1-c/N;
k = N;
while k>0
    if R2_1(k)>0.75
        k = k-1;
    else c = k;
        k = 0;
    end
end
p_3 = 1-c/N;
k = N;
while k>0
    if R2_2(k)>0.75
        k = k-1;
    else c = k;
        k = 0;
    end
end
p_4 = 1-c/N;
t_stat_1 = sort(t_stat_1);
counter1 = 0;
for i = 1:length(t_stat_1) 
    if t_stat_1(i)>1.98 | t_stat_1(i)<-1.98
        counter1 = counter1 +1;
    end 
end
t_stat_2 = sort(t_stat_2);
counter2 = 0;
for i = 1:length(t_stat_2) 
    if t_stat_2(i)>1.98 | t_stat_2(i)<-1.98
        counter2 = counter2 +1;
    end 
end
fprintf('%30s %7.2f %6s %7.2f \n', "The t-test reject the null of nonsignicant impact of one variable on the other for the First scenario", counter1/100,'. For second scenario', counter2/100);
fprintf('%15s %7.4f %6s %7.4f \n', "The relative frequency that the R2 is larger than 0.5 for first Scenario:", p_1,'. For second scenario',p_2);
fprintf('%15s %7.4f %6s %7.4f \n', "The relative frequency that the R2 is larger than 0.75 for first Scenario:", p_3,'. For second scenario',p_4);

%% Please my own functions for R2 and t-stat to make plots 
% R^2
function R2 = calculateR2_mazur(x,y)
% calculateR2_mazur performs regression one variable on the other and
% storing the R^2
% Input:	- x: input x
%           - y: input y
% Output:	- R2: R-squared value 
% USAGE: xx = calculateR2_mazur (x, y)
%
% Author: Aleksander Mazur (SGH), 2022. 

X = [ones(size(x(:,1))),x(:,1)];
Y = y(:,1);
Coeff = X\Y;
yHat = X*Coeff;
res = Y-yHat;
yBar = mean(Y);
regRes = yHat-yBar;
SSR = regRes'*regRes;
SSE = res'*res;
SST = SSR+SSE;
R2 = 1-SSE/SST;
end 

% t-student statistics
function t_stat = calculate_t_stat_mazur(X, Y)
% calculate_t_stat_mazur performs regression one variable on the other and
% calculating the t-statistics 
% Input:	- X: input x
%           - Y: input y
% Output:	- t_stat: t-student statistics 
% USAGE: xx = calculate_t_stat_mazur (X, Y)
%
% Author: Aleksander Mazur (SGH), 2022. 

Coeff = (inv(X'*X))*(X')*Y;
Var = (((Y-X*Coeff)')*(Y-X*Coeff))/(size(Y,1)-size(X,2));
SE = diag((Var.*inv(X'*X))).^(0.5); 
t_stat = Coeff./SE;
end 