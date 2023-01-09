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
