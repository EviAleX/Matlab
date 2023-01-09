% Praca domowa, zadanie 6. Meeting II
clc; % Clear command window 
clear; % Remove items from workspace, freeing up system memory

% Using hint from Lecture2_Ex_5_hint.m.
n=10000;

% Sample from uniform (0,1)
a=rand(n,1);

% Sample from mixture of two Gaussian distributions
dane_symul=(a<0.5).*(randn(n,1)*1+1) + (a>=0.5).*(randn(n,1)*1+4);

%% Subplot #1
% 1. Kernel density 
subplot(2,2,1) % 1 - n. of rows, 3 - n. of colums, 1 - position
[xi,f] = ksdensity_mazur(dane_symul);  
plot(xi, f,'LineWidth', 2,'Color','blue');
hold on

% 2. True density 
mu = [1, 4]; % 1 for first G. distribution, 4 for second G. distribution
sigma = [1, 1]; % -/- 
x = linspace(-10, 10, 500); % 500 stands for number of nsteps, 500 * 2 = 1000 - same length for vectors, ksdensity_mazur

y_1 = 1/((2*pi)^(0.5)*sigma(1))*exp(-((x-mu(1)).^2)/(2*sigma(1)^2)); % density function for first G. distirbution 
y_2 = 1/((2*pi)^(0.5)*sigma(2))*exp(-((x-mu(2)).^2)/(2*sigma(2)^2)); % density function for second G. distirbution 
y = (y_1 + y_2)/2; % whole density 

plot(x, y,'LineWidth', 2,'Color','red')
title({'Nonparametric density', 'estimate (I)'},'FontSize', 14)
legend({'kernel density','true density'})
axis([-5 10 0 0.25]) % function axis allows to set axis limits and aspect ratios: using 4-elements vectors 
hold off

%% Subplot #2
% 1. Histogram density 
subplot(2,2,2) % -/-
histogram(dane_symul, 30, Normalization = "pdf", FaceColor = "blue") % 1 - data, 30 - nbins, to compute density probability we need to use "pdf" normalization
hold on

% 2. True density 
plot(x, y,'LineWidth',2,'Color','red') 
title({'Nonparametric density', 'estimate (II)'},'FontSize', 14)
legend({'histogram density','true density'})
axis([-5 10 0 0.25]) % -/-
hold off

%% Subplot #3
% Frequency histogram for Gaussian PIT
subplot(2,2,3) % -/-    
data_1 = sort(dane_symul); % sorting the data 
PIT_gaussian = normcdf(data_1, mean(data_1), std(data_1)); % calculating PIT assuming standard normal distribution 
h = histogram(PIT_gaussian, 30, Normalization = "probability", FaceColor = "blue");
title({'Frequency histogram', 'for Gaussian PIT'}, 'FontSize', 14)

%% Subplot %4
% Q-Q plot
subplot(2,2,4) % -/-

% set probability levels for quantiles
alpha = 0.0001:0.005:0.9999; 
x_2 = linspace (-3, 3, 1000); 

norm = normalize(dane_symul); % normalizing data 
theoretical_q = norminv(alpha); % Lecture2_Ex3 and Ex4
empirical_q = quantile(norm, alpha); % -/-

scatter(theoretical_q, empirical_q, 'blue', 'filled'); % using scatter plot with circular markers at the locations specified by the theo._q and empr._q
xlabel('standard normal quantiles', 'FontSize', 12);
ylabel('quantiles for std. data', 'FontSize', 12);
hold on
plot(x_2, x_2, 'LineWidth', 3, 'Color', 'red'); 
title('Q-Q plot', 'FontSize', 14);
axis([-2.2 2.2 -3 3])
hold off

%% Please use my own ksdensity function for plot #1 
function [xi, f, h] = ksdensity_mazur (data, nsteps)

% ksdensity_mazur performs kernel density estimation
% Input:	- data: input data, one-dimensional
%           - nsteps: optional number of abscis points. 
% Output:	- xi: equispaced abscis points
%			- f: estimates of p(x)
%           - h: sigma, also called bandwidth 
% USAGE: [xx,pp,dd] = ksdensity_mazur(data, nsteps)
%
% Author: Aleksander Mazur (SGH), 2022. 

if (nargin < 2), nsteps = 1000; end % default number of abscis points
N = length(data); % number of data points 
[minA, maxA] = bounds(data); % calculating the min and max value from the data 
min_max = [minA, maxA]; % storing those values into array 

xi = linspace(min_max(1), min_max(2), nsteps); % linspace stands for linearly spaced vector 
f = zeros(size(xi)); % creating an empty array, which will contain function values 

sigma = std1(data); % stoting std value into variable
h = 1.059*sigma*N^(-1/5); % optimal smoothing parametr for the Gaussian kernel (Silverman , 1986)

% function std1, because I do do my homework
    function y_std = std1(series) % using this funtion to calculate the std 
      y_std = ((ones(1, size(series,1))*((series-mean1(series)).^2))/size(series,1)).^0.5;
    end

% kernel density estimation
for i = 1:N
    K = 1/sqrt(2*pi*h^2)*exp(-(data(i)-xi).^2/(2*h^2)); % calculating the Gaussian kernel
	f = f + 1/N*K; % calculating the kernel density function 
end
end