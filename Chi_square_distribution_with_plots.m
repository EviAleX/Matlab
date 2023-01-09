% Praca domowa, zadanie 20
clc; % Clear command window 
clear; % Remove items from workspace, freeing up system memory

% Szanowna Pani Profesor, w swojej pracy domowej zamieściłem 2 rozwiązania
% dla generacji wykresów "Density functions for chi2 distribution",
% umieściłem w pool "for", proszę przeczytać koniec loop "for": linijki
% 50-53

sz = 100000; % number of values sampled 
matrix_A = [1 3 5 50]; % matrix A with the corresponding degrees of freedom

matrix_B = [15 20 30 100]; % matrix of numbers, which will specify the limits for the current axes
matrix_C = [0.75 0.3 0.27 0.2]; % -/-
matrix_D = [1.3 0.26 0.17 0.05]; % -/-

% Generating figures using "for" loop 
for i=1:4
    % Data for histograms 
    nu = matrix_A(i); % nu - degrees of freedom of corresponding i = 1:4
    r = chi2rnd(nu, sz, 1); % generates an array of random numbers from the chi-square distribution, where sz indicates the size of dimension.
    b = matrix_B(i); % parameters for axis 
    c = matrix_C(i);

    % Data for density functions - Solution 1 
    x1 = 0:0.01:1000; % generation of a sequence of values from 0 to 1000 with a deviation of 0.01. 1000:0.01 = 100000 samples
    y1 = pdf('Chisquare', x1, nu); % probability density function using 'Chisquare' distibution object 
    d = matrix_D(i);

    % Data for density functions - Solution 2
    % x2 = linspace(0, max(r), sz); % generates n points 
    % y2 = 1./((2)^(nu/2).*gamma(nu/2)).*x2.^(nu/2-1).*exp(-x2/2); % chi2 pdf function - was taken from wikipedia 

    % Generating histograms for chi2 distribution with no of degrees of
    % freedom: 1, 3, 5, 50 
    subplot(4,2,i*2-1) % the i-th subplot (corresponding to the i-th replication) is generated; i*2-1 for i=1,2,3,4 -> 1,3,5,7 places where plots will be located 
    h = histogram(r, 22,'Normalization','probability'); % 22 stands for number of bins, 'Normalization' - type of normalization, and 'probability' - relative probability
    xlabel('x');
    title({'Frequency histogram for \chi^2 distribution', ['with no. of degrees of freedom: ' num2str(nu)]}, 'FontSize', 9); % title and corresponding degrees fo freedom
    axis ([0 b 0 c]); % function axis allows to set axis limits and aspect ratios: using 4-elements vectors 
    h.NumBins
    h

    % Generating density functions for chi2 distribution with no of degrees of
    % freedom: 1, 3, 5, 50 - Solution 1 

    subplot(4,2,i*2) % the i-th subplot (corresponding to the i-th replication) is generated; i*2 for i=1,2,3,4 -> 2,4,6,8 places where plots will be located 
    plot(x1, y1, 'LineWidth', 2); % creates a 2D_line plot 
    xlabel('x');
    title({'Density function for \chi^2 distribution', ['with no. of degrees of freedom: ' num2str(nu)]}, 'FontSize', 9); %line 39
    axis ([0 b 0 d]); % line 40

    % Generating density functions for chi2 distribution with no of degrees of
    % freedom: 1, 3, 5, 50 - Solution 2 
    % Proszę zamienić x1 na x2 oraz y1 na y2 w linijce 46 z funkcją plot
    % oraz odblokować linijki 31-32 oraz zablokować 26-27 z pomocą %
end

% Sprawdzić poprawność density function można za pomocą wbudowanej w Matlab
% funckji:
% x_test = 0:0.01:1000;
% y_test = chi2pdf(x_test, 1); % 1 degree of freedom
% figure;
% plot (x_test, y_test, 'LineWidth', 2);
% axis ([0 15 0 1.7]);