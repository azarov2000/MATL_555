function [res] = integrate(f,x0,xend,h)

% f - функция правой части дифференциального уравнения dy/dt = f(t, y)
% tspan - вектор с начальным и конечным временем [t0, tf]
% y0 - начальное значение y(t0)
% h - шаг интегрирования

t = x0:h:xend;
N = length(t);
y0 = 0;
y = zeros(N, length(y0));
y(1, :) = y0;

for i = 1:N-1
    k1 = h * f(t(i), y(i, :));
    k2 = h * f(t(i) + h/2, y(i, :) + k1/2);
    k3 = h * f(t(i) + h/2, y(i, :) + k2/2);
    k4 = h * f(t(i) + h, y(i, :) + k3);

    y(i+1, :) = y(i, :) + (k1 + 2*k2 + 2*k3 + k4)/6;
end
res = y(end)-y(1);





end