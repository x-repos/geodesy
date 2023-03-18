close; clc; clear;
% initial data
h = [0, 5, 10, 15, 20, 25, 30, 35]';
moisture = [124, 78, 54, 35, 30, 21, 22, 18]';

hold
plot(h, moisture,'+')
n = length(moisture);
% use linear leastquare to find coefficients
% of the function of y = a + b*x + c*x^2
A = [ones(n, 1), h, h.^2];
N = A'*A;
coef = N^-1*A'* moisture;
h1 = linspace(min(h), max(h), 100);
f = @(x)coef(1) + coef(2)*x + coef(3)*x.^2;
plot(h1, f(h1))
legend("real data", "regression")
fprintf("coefficents a, b, c are: %5.2f %5.2f %5.2f\n", ...
  coef(1), coef(2), coef(3))
fprintf("prediction at h = 40 is %5.2f\n", f(40))
% variance
meanvalue = sum(h)/n;
variancex = 0;
for i = 1:n
  variancex = variancex + (h(i)-meanvalue)^2;
end
variancex = variancex/n;
g = @(x)coef(2) + 2*coef(3)*x;
variancey = g(meanvalue)^2 * variancex;
fprintf("variance of y %5.2f", variancey)
