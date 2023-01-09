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
