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