function Trapezoidal_Rule
close all; clc;

% initial condition: y(0) = 1
% f: differential equation
% g: solution under initial condition
f = @(x, y) x + y;
g = @(x) 2 * exp(x) - x - 1;

xold = 0;
yold = 1;
h = 0.1;

absdiff = 0;

figure(4);
plot(xold, yold, 'r.');

n = 100;

for i = 1:n
    xnew = xold + h;
    temp = @(ynew) ynew - yold - h/2 * (f(xold, yold) + f(xnew, ynew));
    ynew = fzero(temp, 1);

    drawnow;
    hold on;
    plot(xnew, ynew, 'b.');

    exsol = g(xnew);
    drawnow;
    hold on;
    plot(xnew, exsol, 'k.');

    absdiff = absdiff + abs(ynew - exsol);

    xold = xnew;
    yold = ynew;
end
fprintf("The sum of absolute difference is %d", absdiff);
end
