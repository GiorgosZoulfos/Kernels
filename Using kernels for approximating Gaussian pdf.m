clear; clc; close all;

%% 1000 realizations of U[0,1]

N = 1000;
unif = rand(N,1);

%% Aproximate pdf for different h

h = 0.00001;

M = 3001;
sum = zeros(M,1);

x = -1;
step = 0.001;

for j = 1:M
    
    for i = 1:N
        sum(j) = sum(j) + kernelG(x - unif(i), h);
    end

    sum(j) = sum(j)/N;
    x = x + step;
    
end

plot(-1:.001:2, sum)
title('h = ' + string(h))
grid on


%% Functions used above

% Gaussian kernel
function out = kernelG(x, h)
    
    out = exp(-x^2/(2*h))/sqrt(2*pi*h);
    
end