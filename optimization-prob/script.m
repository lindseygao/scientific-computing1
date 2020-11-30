% Optimization
clear all; close all; clc
% Object attached to spring
% mass: m = 1 kg; spring constant: k = 1/4 N/m; damping force: 1 N
% x(t) function that gives the position of the object at time t where x = 0
% represents the equilibrium position of the spring: x(t) = 1/3*exp(-t/2) +
% 3*t*exp(-t/2)
% x(0) = 1/3; v(0) = 17/6
% a) Finding the maximum value of x(t) and t_max when that occurs using
% fminbnd (0 < t < 10)
x = @ (t) -((1/3)*exp(-t/2) + 3*t*exp(-t/2));
xprime_orig = @ (t) (((-1/6)*exp(-t/2)) + 3*exp(-t/2) - (3/2)*t*exp(-t/2));
xpp_orig = @ (t) (1/12)*exp(-t/2) - 3*exp(-t/2) + (3/4)*t*exp(-t/2)
t_max = fminbnd(x,0,10);
ans1 = [t_max -x(t_max)];
% b) Use golden section search to find t_max & stop when b - a < 1e-3
a = 0;
b = 10;
c = (-1 + sqrt(5)) / 2;
tol = 1e-3;
t1 = (b-a)*(1-c);
t2 = (b-a)*c;
iterations = 0;
while abs(b - a) > tol && iterations < 100
    if x(t1) < x(t2)
        b = t2;
        t2 = t1;
        t1 = a + (b-a)*(1-c);
    elseif x(t1) > x(t2)
        a = t1;
        t1 = t2;
        t2 = a + (b-a)*c;
    else
        a = t1;
        b = t2;
    end
    iterations = iterations + 1;
end
ans2 = [a b];
% c) Use Newton's method to find t_max by finding the root of x'(t)

[ans3,ans4] = Newton(xprime_orig,xpp_orig,0,1e-3,100);


%% Problem 2
f = @ (x,y) (2-x).^2 + (y - x.^2).^2;
f = @ (p) f(p(1),p(2)); %
p = [0; 5]; % initial guess
ans5 = fminsearch(f, p);
% Finding the minimum using gradient descent
fgrad = @ (x,y) [-2*(2-x) - 4*x*(y-x.^2);
                2*(y-x.^2)];
fgrad = @ (p) fgrad(p(1),p(2)); %
iter = 0;
tol = 1e-4;
grad = fgrad(p);
tic
while norm(grad,Inf) > tol && iter < 10000
    grad = fgrad(p);
    phi = @ (t) p - t*grad;
    f_of_phi = @ (t) f(phi(t));
    tmin = fminbnd(f_of_phi,0,1);
    p = phi(tmin);
    iter = iter + 1;
end
ans6 = p;
grad_min = toc
ans7 = iter - 1;

% Using a fixed step size for gradient descent
maxIter = 10000;
iter = 0;
p = [0;5];
grad = fgrad(p);
tstep = .01;
tic
while norm(grad,Inf) > tol && iter < maxIter
        grad = fgrad(p);
        p = p - tstep*grad;
        iter = iter + 1;
end
ans8 = iter-1;
grad_step = toc
iter = 0;
tstep = .03;
p = [0;5];
grad = fgrad(p);
while norm(grad,Inf) > tol && iter < maxIter
        grad = fgrad(p);
        p = p - tstep*grad;
        iter = iter + 1;

end
ans9 = iter-1;
% step size change to .05
iter = 0;
tstep = .05;
p = [0;5];
grad = fgrad(p);
while norm(grad,Inf) > tol && iter < maxIter
        grad = fgrad(p);
        p = p - tstep*grad;
        iter = iter + 1;
end
ans10 = iter-1;
%%
function [xk,iteration] = Newton(f,fprime,xk,tol,maxIter)
% Performs Newton's method
% Inputs:
%   f = function to find the roots of
%   fprime = derivative of f
%   xk = initial guess
%   tol = tolerance for stopping criteria
%   maxIter = maximum number of iterations

change = 2*tol;
iteration = 0;
while change > tol && iteration < maxIter
    xkplus1 = xk - f(xk)/fprime(xk); 
    change = abs(xkplus1-xk);
    xk = xkplus1; 
    iteration = iteration+1;
end

end  