% Initial condition: y(1) = 3
f = @(x, y) sin(x)/x;
g = @(x) sinint(x) - sinint(1) + 3;

% Initial condition: y(1) = 3
% f = @(x, y) sin(x) + cos(x);
% g = @(x) sin(x) - cos(x) + 3 - sin(1) + cos(1);


prompt = "Please input the method by typing the number: \n " + ...
    "1. Forward_Euler, " + ...
    "\n 2. Backward_Euler, " + ...
    "\n 3. Runge_Kutta, " + ...
    "\n 4. Trapezoidal," + ...
    "\n 5. Adams Bashforth. \n";
method = input(prompt);

% log_2(delta_t), log_2(e(delta_t))
delta_t_list = [1, 0.5, 0.25, 0.125, 0.0625]; 
% delta_final is used to record the max error when the delta_t is 1, 0.5, 0.25, 0.125, 0.0625
error_final = zeros(5);

switch method
   case 1
       for i = 0:4
           [absdiff_FE, e, delta_t] = Forward_Euler(delta_t_list(i+1), 10*(2^i), f, g);
           error_final(i + 1) = e;
       end

   case 2
       for i = 0:4
           [absdiff_BE, e, delta_t] = Backward_Euler(delta_t_list(i+1), 10*(2^i), f, g);
           error_final(i + 1) = e;
       end

   case 3
       for i = 0:4
           [absdiff_RK4, e, delta_t]= Runge_Kutta_4(delta_t_list(i+1), 10*(2^i), f, g);
           error_final(i + 1) = e;
       end

   case 4
       for i = 0:4
           [absdiff_RK4, e, delta_t]= Trapezoidal_Rule(delta_t_list(i+1), 10*(2^i), f, g);
           error_final(i + 1) = e;
       end

   case 5
       for i = 0:4
           [absdiff_RK4, e, delta_t]= Adams_Bashforth_2(delta_t_list(i+1), 10*(2^i), f, g);
           error_final(i + 1) = e;
       end

   otherwise
      fprintf("Invalid input number")
end 

prompt2 = "\n Convergence plot? y/n: ";

conver = input(prompt2, "s");

switch conver
    case 'y'
        convergence_plot(error_final, delta_t_list);
    otherwise
        fprintf("No plot request");
end

% Forward_Euler
function [absdiff,error_max,delta_t] = Forward_Euler(h, n, f, g)
close all; clc;

xold = 1;
yold = 3;

absdiff = 0;

figure(1);
plot(xold, yold, 'r.');

delta_t = 0;
error_max = 0;
for i = 1:n
    xnew = xold + h;
    ynew = yold + h * f(xold, yold);
   
    drawnow;
    hold on;
    plot(xnew, ynew, 'b.');

    exsol = g(xnew);
    drawnow;
    hold on;
    plot(xnew, exsol, 'k.');

    
    drawnow;
    hold on;
    plot(xnew, abs(ynew - exsol));

    absdiff = absdiff + abs(ynew - exsol);
    
    delta_t = h;
    
    error = abs(ynew - exsol);

    if error > error_max
        error_max = error;
    end

    xold = xnew;
    yold = ynew;
end
% fprintf("The sum of absolute difference is %d", absdiff);
end

% Runge_Kutta_4
function [absdiff,error_max,delta_t] = Runge_Kutta_4(h, n, f, g)
close all; clc;


xold = 1;
yold = 3;

absdiff = 0;

figure(3);
plot(xold, yold, 'r.');

delta_t = 0;
error_max = 0;

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

    absdiff = absdiff + abs(ynew - exsol);

    delta_t = h;
    
    error = abs(ynew - exsol);

    if error > error_max
        error_max = error;
    end

    xold = xnew;
    yold = ynew;
end
% fprintf("The sum of absolute difference is %d", absdiff);
end

% Backward_Euler
function [absdiff,error_max,delta_t] = Backward_Euler(h, n, f, g)
close all; clc;

xold = 1;
yold = 3;

absdiff = 0;

error_max = 0;
delta_t = h;

% figure(2);
% plot(xold, yold, 'r.');

for i = 1:n
    xnew = xold + h;
    temp = @(ynew) ynew - yold - h * f(xnew, ynew);
    ynew = fzero(temp, 1);

    drawnow;
    hold on;
    plot(xnew, ynew, 'b.');

    exsol = g(xnew);
    drawnow;
    hold on;
    plot(xnew, exsol, 'k.');

    absdiff = absdiff + abs(ynew - exsol);
    
    error = abs(ynew - exsol);

    if error > error_max
        error_max = error;
    end

    xold = xnew;
    yold = ynew;
end
fprintf("The sum of absolute difference is %d", absdiff);
end

% Trapezoidal_Rule
function [absdiff,error_max,delta_t] = Trapezoidal_Rule(h, n, f, g)
close all; clc;

xold = 1;
yold = 3;

absdiff = 0;

error_max = 0;
delta_t = h;

% figure(4);
% plot(xold, yold, 'r.');

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

    delta_t = h;
    
    error = abs(ynew - exsol);

    if error > error_max
        error_max = error;
    end

    xold = xnew;
    yold = ynew;
end
% fprintf("The sum of absolute difference is %d", absdiff);
end

% Adams_Bashforth_2
function [absdiff,error_max,delta_t] = Adams_Bashforth_2(h, n, f, g)
close all; clc;

xolder = 1;
yolder = 3;

xold = xolder + h;
yold = yolder + h * f(xolder, yolder);

absdiff = abs(yold - g(xold));

% figure(5);
% plot(xolder, yolder, 'r.');
% hold on;
% plot(xold, yold, 'b.');
% hold on;
% plot(xold, g(xold), 'k.');

delta_t = h;
error_max = 0;
for i = 1:(n-1)
    xnew = xold + h;
    ynew = yold + 3*h/2 * f(xold, yold) - h/2 * f(xolder, yolder);

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
    
    delta_t = h;
    
    error = abs(ynew - exsol);

    if error > error_max
        error_max = error;
    end

    xold = xnew;
    yold = ynew;
end
% fprintf("The sum of absolute difference is %d", absdiff);
end

function convergence_plot(e, delta_t)
% figure()
tiledlayout(1,2)

% sp(1) = subplot(1, 1, 1);
nexttile;
plot(delta_t, e);
title('Convergence Plot')

% sp(2) = subplot(2, 2, 2);

log_delta = log2(delta_t);

log_e = log2(e);

nexttile;

plot(log_delta, log_e);
title('Rate of Convergence')

slope = (log_e(1) - log_e(5)) / (log_delta(1) - log_delta(5));

fprintf("The slope is  %d", slope);
end