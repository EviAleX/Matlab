% Praca domowa, zadanie 5. Meeting II
clc; % Clear command window 
clear; % Remove items from workspace, freeing up system memory

% Using hint from Lecture2_Ex_5_hint.m.
n=10000;

% Sample from uniform (0,1)
a=rand(n,1);

% Sample from mixture of two Gaussian distributions
dane_symul=(a<0.5).*(randn(n,1)*1+1) + (a>=0.5).*(randn(n,1)*1+4);

% Making variables for True density
mu = [1 4];
sigma = [1 1]; 
x = -10:0.01:10;

% True density itself
y_1 = 1/((2*pi)^(0.5)*sigma(1))*exp(-((x-mu(1)).^2)/(2*sigma(1)^2)); % true density for first G. distribution
y_2 = 1/((2*pi)^(0.5)*sigma(2))*exp(-((x-mu(2)).^2)/(2*sigma(2)^2)); % true density for second G. distribution
y = (y_1 + y_2)/2; % whole density

% Making data for Kernel density
N = length(dane_symul); % number of data points 
[minA, maxA] = bounds(dane_symul); % calculating the min and max value from the data 
min_max = [minA, maxA]; % storing those values into array 

xi = linspace(min_max(1), min_max(2), n); % linspace stands for linearly spaced vector 
f = zeros(size(xi)); % creating an empty array, which will contain function values 
sigma = std1(dane_symul); % stoting std value into variable

% Smoothing parametrs h 
h_1 = 1.059*sigma*N^(-1/5); % optimal smoothing parametr for the Gaussian kernel (Silverman , 1986)
h_2 = h_1*5; % optimal smoothing parametr for middle figure
h_3 = 0.15*h_1; % optimal smoothing parametr for right figure 
h_box = [h_1 h_2 h_3]; % matrix, where stored all h. 

% Making loop to not repeat yourself
for j = 1:3 % j used for smoothing parametr 
    h = h_box(j);
    f = zeros(size(xi));
    K = zeros(1, n);
    for i = 1:N
         K = 1/sqrt(2*pi*h^2)*exp(-(dane_symul(i)-xi).^2/(2*h^2)); % calculating the Gaussian kernel
         f = f + 1/N*K; % calculating the kernel density function 
    end
    subplot(1,3,j);
    plot(xi, f,'LineWidth', 1,'Color','blue'); % kernel density 
    hold on
    plot(x, y,'LineWidth', 1,'Color','red'); % true density 
    legend({'kernel density','true density'}, 'Fontsize', 6); 
    axis([-5 10 0 0.25]) % function axis allows to set axis limits and aspect ratios: using 4-elements vectors 
    set(gcf, 'Position', [100 100 1400 400]) 
    hold off
end 