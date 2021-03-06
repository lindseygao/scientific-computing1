clear all; close all; clc
% Population Data
load('SeaPopData.mat'); % Seattle population data from 
% https://en.wikipedia.org/wiki/Demographics_of_Seattle
% & https://www.seattle.gov/opcd/population-and-demographics

% Finding the line of best fit: P = mt + b where t is years since 1860
lin_fit = polyfit(t,Seattle_Pop,1); % slope of line
lin_est = polyval(lin_fit,2019-1860); % population estimate for 2019

% Finding the best fit quadratic function
quad_fit = polyfit(t,Seattle_Pop,2); % coefficients of quadratic
quad_est = polyval(quad_fit,2019-1860); % 2019 population estimate

% Finding best fit degree 5 polynomial
deg5 = polyfit(t,Seattle_Pop,5);
deg5_est = polyval(deg5,2019-1860); % 2019 pop estimate

% Finding best fit degree 9 polynomial
deg9 = polyfit(t,Seattle_Pop,9);
deg9_est = polyval(deg9,2019-1860); % 2019 pop estimate

% Plotting Seattle population data & all best fit polynomials from above
plot(t,Seattle_Pop,'ko', 'Linewidth', [1]), hold on
t_plot = 0:.01:160;
lin_plot = polyval(lin_fit,t_plot);
plot(t_plot,lin_plot,'b', 'Linewidth', [1])
quad_plot = polyval(quad_fit,t_plot);
plot(t_plot,quad_plot,'r', 'Linewidth', [1])
deg5_plot = polyval(deg5,t_plot);
plot(t_plot,deg5_plot,'m', 'Linewidth', [1])
deg9_plot = polyval(deg9,t_plot);
plot(t_plot,deg9_plot,'g', 'Linewidth', [1])
xlabel('Years since 1860')
ylabel('Seattle Population')
legend('data', 'deg 1', 'deg 2', 'deg 5', 'deg 9', 'Location','Best','Fontsize', [10])
ylim([0 800000])
print -dpng Problem2.png

%% Atmospheric CO2 Data measured at the Mauna Loa observatory in Hawaii
% Calculating exponential fit using least-squares fitting method
load('CO2_data.mat') % monthly CO2 averages since 1958
% y = a*exp(rt)
% ln(y) = ln(aexp(rt)) = ln(a) + rt
ln_lin = polyfit(t,log(CO2),1)
ans5 = ln_lin(1) % r in rt (slope)
ans6 = exp(ln_lin(2)) % y-int

% Potting the CO2 data & previous exponential fit
log_lin = ans6*exp(ans5*t);
t_plot = 0:.01:65
plot(t,CO2, '-k.'), hold on
plot(t_plot, exp(polyval(ln_lin,t_plot)), 'r', 'Linewidth', 2)
xlabel('Years since 1958')
ylabel('Atmospheric CO_2')
title('Data Linearization Method for Exponential Fit')
xlim([0 65])
legend('data', 'fit curve', 'Location', 'Best')
print problem4.png -dpng

%% Altering the exponential fit by shifting it up by a constant
% y = a*exp(rt) => y = a*exp(rt) + b
% Does this by minimizing the sum of squared errors (root-mean squared
% error) & using fminsearch()
load('CO2_data.mat');
% initial guess: a = 30, r = .03, b = 300
adapter = @ (p) sumSquaresError(p(1),p(2),p(3));
p0 = [30; .03; 300];
[ans7 ans8] = fminsearch(adapter, p0) % ans7 = optimal parameters: [a; r; b]
% ans8 = min value of sum of squared errors


% Plot that contains the atmospheric CO2 data & best fit curve
t_plot = 0:.01:65;
plot(t,CO2, '-k.'), hold on
plot(t_plot, ans7(1)*exp(ans7(2)*t_plot)+ans7(3), 'r', 'Linewidth', 2)
xlabel('Years since 1958')
ylabel('Atmospheric CO_2')
xlim([0 65])
legend('data', 'fit curve', 'Location', 'Best')
title('Best Fit Exponential Curve Minimizing Sum of Squared Errors')
print problem6.png -dpng

%% Altering the best fit curve to take into consideration oscillations
% y = a*exp(rt) + b => y = a*exp(rt) + b + c*sin(d(t-e))
% Still minimizing sum of squared errors
adapter1 = @ (g) SSEoscillation(g(1),g(2),g(3),g(4),g(5),g(6));
g0 = [ans7(1); ans7(2); ans7(3); 3; 6; 0];
[ans9 ans10] = fminsearch(adapter1, g0)
% Plotting best fit curve with oscillations against CO2 data
t_plot = 0:.01:65;
plot(t,CO2, '-k.'), hold on
plot(t_plot, ans9(1)*exp(ans9(2)*t_plot)+ans9(3)+ans9(4)*sin(ans9(5)*(t_plot-ans9(6))), 'r', 'Linewidth', 1)
xlabel('Years since 1958')
ylabel('Atmospheric CO_2')
xlim([0 65])
legend('data', 'fit curve', 'Location', 'Best')
title('Best Fit Oscillation Curve Minimizing Sum of Squared Errors')
print problem8.png -dpng

% function for oscillations
function [error] = SSEoscillation(a,r,b,c,d,e)
y = @(t) a*exp(r*t)+b+c*sin(d*(t-e));
load('CO2_data.mat');
squared_errors = (y(t)-CO2).^2;
error = sum(squared_errors);
end

% calculates sum of squared errors
function [error] = sumSquaresError(a,r,b)
y = @(t) a*exp(r*t)+b;
load('CO2_data.mat');
squared_errors = (y(t)-CO2).^2;
error = sum(squared_errors); % connected with the 2-norm
end