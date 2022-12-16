function Adams_Bashforth_2
close all; clc;

% initial condition: y(0) = 1
% f: differential equation
% g: solution under initial condition
f = @(x, y) x + y;
g = @(x) 2 * exp(x) - x - 1;

xolder = 0;
yolder = 1;
h = 0.1;

xold = xolder + h;
yold = yolder + h * f(xolder, yolder);

absdiff = abs(yold - g(xold));

figure(5);
plot(xolder, yolder, 'r.');
hold on;
plot(xold, yold, 'b.');
hold on;
plot(xold, g(xold), 'k.');

n = 100;

for i = 1:(n-1)
    xnew = xold + h;
    temp = @(ynew) ynew - yold - 3*h/2 * f(xold, yold) + h/2 * f(xolder, yolder);
    ynew = fzero(temp, 0);

    drawnow;
    hold on;
    plot(xnew, ynew, 'b.');

    exsol = g(xnew);
    drawnow;
    hold on;
    plot(xnew, exsol, 'k.');

    absdiff = absdiff + abs(ynew - exsol);

    xolder = xold;
    yolder = yold;

    xold = xnew;
    yold = ynew;
end
fprintf("The sum of absolute difference is %d", absdiff);
end
