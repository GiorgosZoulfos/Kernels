%% Problem 3.2
clear; clc; close all;

%% Load data from data32.mat

load('data32.mat')
N = length(circles);

lamda = 10;

%% Define the needed matrices

% Labels of data: 1 for stars, -1 for circles
labels = [ones(N,1); -ones(N,1);];

% Matix K: Ki,j = K(xi,xj)
data = [stars; circles];
K = zeros(2*N, 2*N);

for i = 1:2*N
    for  j = 1:2*N
        K(i,j) = kernelG(data(i,:), data(j,:));
    end
end

%% Solve the linear system

A = (K'*K) + lamda*K;
b = K*labels;

% vector coef contains the desired coeficients
coef = linsolve(A,b);

%% Calculate Error

errCircles = 0;
errStars = 0;

for i = 1:N
    
    if  phi(stars(i,:), coef(1:N), coef(N+1:end), N, stars, circles) < 0
        errStars = errStars + 1;
    end
    
    if phi(circles(i,:), coef(1:N), coef(N+1:end), N, stars, circles) > 0
        errCircles = errCircles+ 1;
    end
    
end

errStars/N
errCircles/N
(errCircles + errStars)/(2*N)


%% Compute border in 2-D space

allX = -1.05:0.001:1.5;
allY = -.2:0.001:1.2;

points = zeros(length(allX), 2);

for i = 1:length(allX)
    
    min = 10^9;
    ind = 0;
    
    for j = 1:length(allY)
        phix = phi([allX(i), allY(j)], coef(1:N), coef(N+1:end), N, stars, circles);
        if phix < min
            ind = j;
            min = phix;
        end
    end
   
    points(i,1) = allX(i);
    points(i,2) = allY(ind);
    
end

%% Plot data

labels = [ones(21,1); -ones(21,1)];
gscatter([stars(:,1); circles(:,1)],[stars(:,2); circles(:,2)], labels)
hold
plot(points(:,1), points(:,2))
title('L = ' +string(lamda) + ', h = 0.0001')

%% Functios used above

% Error function
function out = errorFun(in, N, lamda, stars, circles, A, B)

    out = lamda * norm(phi(in, A, B, N, stars, circles))^2;
    for i = 1:N
        out = out + (1 - phi(stars(i,:), A, B, N, stars, circles))^2 ...
            + (1 + phi(circles(i,:), A, B, N, stars, circles))^2;
    end

end


% Derivative of error function
function out = errorFunDer(in, N, lamda, stars, circles, A, B)

    out = 2*lamda*phi(in, A, B, N, stars, circles).*phiDer(in, N, stars, circles);
    
    for i = 1:N
        out = out - 2*(1 - phi(stars(i,:), A, B, N, stars, circles)).*phiDer(stars(i,:), N, stars, circles) ...
            + 2*(1 + phi(circles(i,:), A, B, N, stars, circles)).*phiDer(circles(i,:), N, stars, circles);
    end

end

% Phi function 
function out = phi(in, A, B, N, stars, circles)

    out = 0;
    
    for i = 1:N
        out = out + A(i)*kernelG(in, stars(i,:)) + B(i)*kernelG(in, circles(i,:));
    end

end

% Derivative of phi 
function out = phiDer(in, N, stars, circles)

    out = zeros(2*N,1);
    for i = 1:N
        out(i) = kernelG(in, stars(i,:));
        out(i+N) = kernelG(in, circles(i,:));
    end
    
end

% Gaussian Kernel
function out = kernelG(in, stCi)

    h = 0.0001;
    out = exp((-1/h)*(norm(in - stCi)^2));

end


    


