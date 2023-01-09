% Praca domowa, zadanie 2
clc; % Clear command window 
clear; % Remove items from workspace, freeing up system memory

% Creating normal CDF
daily = xlsread('Apple_daily_prices.xlsx'); % Reading the first worksheet from Excel
logreturns_daily = log(daily(2:end,5))-log(daily(1:end-1,5)); % Matrix of financial log returns can be easily defined as

n = length(logreturns_daily); % Calculating the lenght of vector 

x = min(logreturns_daily):(abs(min(logreturns_daily))+abs(max(logreturns_daily)))/(n-1):max(logreturns_daily); % Making a sorted vector from min value to max value of x, to calculate normcdf
y_normcdf = normcdf (x, mean(logreturns_daily), std(logreturns_daily)); % creating y for normcdf function 

returns = sort(logreturns_daily); % Sorting the vector 

% Creating empirical CDF
y_cdf = zeros(n,1); % Creating an empty matrix to fill up with values 
k = 1;

for i = 1:n % Making the vector to create empirical cdf, the same lenght as normcdf
    while returns(k) < x(i)
        k = k + 1;
    end
    y_cdf(i) = (k-1)/n;
end

% Making plots
plot(x, y_normcdf, 'LineWidth', 2, 'Color', 'blue') % Making a plot for normcdf
hold on; % Telling plot to stop
plot(x, y_cdf,'LineWidth',2,'Color','red') % Making a plot for empirical cdf
set(gca, 'FontSize', 14); % Setting the format size
xlabel('return');
ylabel('CDF');
legend({'Normal CDF','Empirical CDF'}, 'Location', 'northwest'); % Making a legend
set(gca,'xTick',-0.1:0.05:0.15) % Setting x-axis tick values 
axis([-0.15 0.15 0 1]) % function axis allows to set axis limits and aspect ratios: using 4-elements vectors
legend boxoff 
hold off; 

% Calculating Lilliefors test
ks = max(y_cdf'- y_normcdf); % sup = max 
g = (0.83+n)/sqrt(n)-0.01; % Calculating the g(n) value, page 31 from script 
p = 0.895/g; % Calculating the p-value for alfa = 0.05

fprintf('%12s %7.4f \n', "Lilliefors test statistic is:", ks); % Using well known fprintf fuction
alfa = 0.05;
if (p < alfa) % Basic if/else statement 
    fprintf('%15s %6.4f %13s \n', "Critical value is", p,", thus we reject the H_0.");
else
    fprintf('%15s %6.4f %15s \n', "Critical value is", p,", thus fail to reject the H_0.");
end

% Additionaly to compare results we can use Matlab build-in function
% [h,p,kstat,critval] = lillietest( );