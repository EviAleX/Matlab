% Praca domowa, zadanie 1
clc; % Clear command window 
clear; % Remove items from workspace, freeing up system memory

daily = xlsread('Apple_daily_prices.xlsx'); % Reading the first worksheet from Excel 
weekly = xlsread('Apple_weekly_prices.xlsx');
monthly = xlsread('Apple_monthly_prices.xlsx');

logreturns_daily = log(daily(2:end,5))-log(daily(1:end-1,5)); % Matrix of financial log returns can be easily defined as
logreturns_weekly = log(weekly(2:end,5))-log(weekly(1:end-1,5));
logreturns_monthly = log(monthly(2:end,5))-log(monthly(1:end-1,5));

statistics = zeros(6,3); % Creating an empty matrix to fill up with values 
matrix_for_statistics = {logreturns_daily logreturns_weekly logreturns_monthly}; % Specified cell matrix filled by created alredy variables, easier to use for loop

for i=1:3 % Basic loop filled with statistics, which then will be filled in matrix created earlier 
    k = matrix_for_statistics{i}; 
    statistics(1,i) = mean(k); % Calculating secriptive statistics according to the task 
    statistics(2,i) = std(k);
    statistics(3,i) = skewness(k);
    statistics(4,i) = kurtosis(k);
    [statistics(5,i), statistics(6,i)] = jbtest1(k);
end

indices = {'daily return', 'weekly returns', 'monthly returns'}; % Making an array with titles 

% Display the outcome 
disp('--------------------------------------------------------------------------------------');
disp('                   Mean   Std.Dev.  Skewness  Kurtosis    J-B stat      J-B p-value'   );
disp('--------------------------------------------------------------------------------------');
for i=1:3 % Making a loop to fill up statistics into the display outcome, copied from example of page 24
    fprintf('%15s %8.3f %8.3f %8.3f %10.3f %13.3f %11.3f \n', indices{i}, statistics(:,i));
end 

% Nie ma statystycznych podstaw powiedzić że te rozkłady pochodzą z rozkłądu normalnego. 
% Wielka kurtoza wskazuje że w tych rozkłądach możemy się spodziewać grubych ogonów. 
% Dodatkowo wysoka dodatnia wartość dla kurtozy wskazuje na to, iż krańce rozkładu są dłuższe niż te dla rozkładu normalnego
% Ujemna skośność mówi o gwałtowniejszych spadkach w tej akcji.