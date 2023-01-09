% Praca domowa, zadanie 2. Meeting V
clc; % Clear command window 
clear; % Remove items from workspace, freeing up system memory

n = 1000; % n. of observations
matrix_a = [0.2 0.7 0.05 0.01 0.01]; % matrix of parameters alpha
matrix_b = [0.7 0.2 0.8 0.95 0.98]; % matrix of parameters beta
w = zeros(5, 1); % empty matrix for w 
for i = 1:5 
    w(i) = 0.3/(1-(matrix_a(i)+matrix_b(i))); % E(at^2)=w/(1-(alpha+beta)) - GARCH(1,1)
end

% creating a loop 
for i=1:5
    rng(11)
    a = matrix_a(i); % each parameter alpha for each iteration
    b = matrix_b(i); % each parameter beta for each iteration
    sigma_square = zeros(n,1); % empty matrix for conditional variances
    a_simple = zeros(n,1); % empty matrix for returns
    eps = randn(n,1); 
    sigma_square(1,1) = 0.3; % to initiliaze loop
    a_simple(1,1) = 1; % -/-
    for j=1:n % according to the task 
        a_simple(j) = eps(j)*(sigma_square(j)^0.5);
        sigma_square(j+1) = w(i) + a*(a_simple(j)^2) + b*sigma_square(j);
    end 
    subplot(5,2,i*2-1) % creating plots
    plot(sigma_square(1:end-1), 'LineWidth', 2, 'Color', 'r');
    title(['\alpha=' num2str(a) ', \beta=' num2str(b)], 'FontSize', 10);
    subplot(5,2,i*2) % creating plots
    plot(a_simple);
    title(['\alpha=' num2str(a) ', \beta=' num2str(b)], 'FontSize', 10);
end
