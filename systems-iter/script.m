clear all; close all; clc

% Implementing discritzed version of Poisson's equation: -x"(t) = b(t)
% a) constructing 80x80 matrix (2 on diagonal & -1 on diagonals above and
% below main diagonal)
% b(j) = exp(-15pi/10)*sin(15*pi*j/81)
diagonal(1:80, 1) = 2;
A80_diag = diag(diagonal);
diag1(1:79, 1) = -1;
A80 = A80_diag + diag(diag1, 1) + diag(diag1, -1);
j = 1:80;
b(j) = exp(-15*pi/10)*sin((15*pi*j)/81);
b = b.';
% b) exact solution to Ax = b
ans1 = A80\b; % exact solution
% c) Implementing jacobi method to above system
D = diag(diag(A80));
T = A80 - D;
Mj = -D\T;
c = D\b;
lambda = eig(Mj);
ans2 = max(abs(lambda));
% d) using jacobi method to solve Ax = b
int_guess = diag(A80);
tol = 1e-4;
maxIter = 10000;
[x_jacobi, iter_jacobi] = Jacobi(A80, b, int_guess, tol, maxIter);
error_jacobi = norm(ans1-x_jacobi, 2);
ans3 = [iter_jacobi error_jacobi];
% e) implementing gauss-seidel method
S = tril(A80);
T = A80 - S;
Mg = -S\T;
c = S\b;
lambda_g = eig(Mg);
ans4 = max(abs(lambda_g));
% f) using gauss-seidel method to solve Ax= b
[x_gauss, iter_gauss] = GaussSeidel(A80,b,int_guess,tol,maxIter);
error_gauss = norm(ans1-x_gauss, 2);
ans5 = [iter_gauss error_gauss];

% Implementing successive over-relaxation method to above system
% a)
w = 1.2;
L = S - D;
U = triu(A80) - D;
P = (1/w)*D + L;
T = ((w-1)/w)*D + U;
M = -P\T;
ans6 = norm(eig(M), inf);
% b) using SOR to solve for x 
[x_SOR, iter_SOR] = SOR(A80, b, int_guess,tol,1.2,maxIter);
error_SOR = norm(ans1-x_SOR, 2);
ans7 = [iter_SOR error_SOR];
% c) Finding the optimal value of w (between 1 & 2) so SOR converges faster
% Optimal w when eigenvalue of M has the largest absolute value (so use the
% smallest of these eigenvalues in absolute value)
currentMax_eig = 1000;
for w = 1:.01:2;
    P = (1/w)*D + L;
    T = ((w-1)/w)*D + U;
    M = -P\T;
    if max(abs(eig(M))) < currentMax_eig
        currentMax_eig = max(abs(eig(M)));
        opt_w = w;
    end
end

ans8 = opt_w; % optimal w
ans9 = currentMax_eig;
% d) Use SOR to solve for x with optimal w
[x_opt, iter_opt] = SOR(A80, b, int_guess,tol,ans8,maxIter);
error_opt = norm(ans1-x_opt, 2);
ans10 = [iter_opt error_opt];

% Solving Ax = b again but calculating the error each time and plotting the
% error vector on log(y) scale against number of iterations vector
% k = num iterations vector = (0 1 2 ...)
% Calculating eigenvalue of M that has the largest in absolute value and
% plotting this (lambda^0 lambda^1 lambda^2 ...) - lambda^k
D = diag(diag(A80));
T = A80-D;
M = -D\T;
c = D\b;
change = 2*tol;
iterations = 0;
xk = int_guess;
lambda_max(1) = max(abs(eig(M)))^iterations;
error_vec(1) = norm(ans1-xk,2);
while change > tol && iterations < maxIter
    xkplus1 = M*xk + c;
    change = norm(xkplus1-xk,Inf); % largest value in absvalue 
    iterations = iterations + 1;
    xk = xkplus1;
    error_vec(iterations + 1) = norm(ans1-xk,2);
    lambda_max(iterations + 1) = max(abs(eig(M)))^(iterations + 1);
end 

k = 0:iterations;
semilogy(k, error_vec, 'r', 'Linewidth', [2]), hold on
semilogy(k, lambda_max, 'k--', 'Linewidth', [2])
xlabel('Iterations, k', 'Fontsize',[12])
ylabel('Error (log scale)', 'Fontsize', [12])
title(['Jacobi Method Error'], 'Fontsize', [15])
legend('Jacobi Error', '\lambda_{max}^{k}', 'Location', 'Best', 'Fontsize', [10])
print -dpng problem4.png

function [x, numIter] = SOR(A,b,x0,tol,w,maxIter)
D = diag(diag(A));
S = tril(A);
L = S - D;
U = triu(A) - D;
P = (1/w)*D + L;
T = ((w-1)/w)*D + U;
M = -P\T;
c = P\b;
change = 2*tol;
iterations = 0;
xk = x0;
while change > tol && iterations < maxIter
    xkplus1 = M*xk + c;
    change = norm(xkplus1-xk,Inf); % largest value in absvalue 
    iterations = iterations + 1;
    xk = xkplus1;
end 
x = xk;
numIter = iterations;
end


function [x,numIter] = GaussSeidel(A,b,x0,tol,maxIter)
S = tril(A);
T = A-S;
M = -S\T;
c = S\b;
change = 2*tol;
iterations = 0;
xk = x0;
while change > tol && iterations < maxIter
    xkplus1 = M*xk + c;
    change = norm(xkplus1-xk,Inf); % largest value in absvalue 
    iterations = iterations + 1;
    xk = xkplus1;
end 
x = xk;
numIter = iterations;
end
    

function [x,numIter] = Jacobi(A,b,x0,tol,maxIter)
D = diag(diag(A));
T = A-D;
M = -D\T;
c = D\b;
change = 2*tol;
iterations = 0;
xk = x0;
while change > tol && iterations < maxIter
    xkplus1 = M*xk + c;
    change = norm(xkplus1-xk,Inf); % largest value in absvalue 
    iterations = iterations + 1;
    xk = xkplus1;
end 
x = xk;
numIter = iterations;
end





