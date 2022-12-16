function absdiff = Runge_Kutta_4(f, g)
close all; clc;

% initial condition: y(0) = 1
% f: differential equation
% g: solution under initial condition
% f = @(x, y) x + y;
% g = @(x) 2 * exp(x) - x - 1;

xold = 0;
yold = 1;
h = 0.1;

absdiff = 0;

figure(3);
plot(xold, yold, 'r.');

n = 100;

for i = 1:n
    k1 = f(xold, yold);
    k2 = f(xold + h/2, yold + h/2 * k1);
    k3 = f(xold + h/2, yold + h/2 * k2);
    k4 = f(xold + h, yold + h * k3);

    xnew = xold + h;
    ynew = yold + (h/6) * (k1 + 2*k2 + 2*k3 + k4);

    drawnow;
    hold on;
    plot(xnew, ynew, 'b.');

    exsol = g(xnew);
    drawnow;
    hold on;
    plot(xnew, exsol, 'k.');

    stem(i, abs(ynew - exsol));

    absdiff = absdiff + abs(ynew - exsol);

    xold = xnew;
    yold = ynew;
end
% fprintf("The sum of absolute difference is %d", absdiff);
end