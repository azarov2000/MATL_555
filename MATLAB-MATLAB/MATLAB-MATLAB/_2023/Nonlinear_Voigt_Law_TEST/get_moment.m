function [interp_t,Sum_m] = get_moment(t,ksi_x,ksi_y,ksi_x_dot,ksi_y_dot,m,N)

% Интерполяция по точкам с равномерным шагом
interp_t = linspace(0,t(end),length(t));
for j=1:length(ksi_x(1,:))
    parfor i=1:length(interp_t)
        x_inrp(i,j) = interp1(t,ksi_x(:,j),interp_t(i));
        y_inrp(i,j) = interp1(t,ksi_y(:,j),interp_t(i));
        x_dot_inrp(i,j) = interp1(t,ksi_x_dot(:,j),interp_t(i));
        y_dot_inrp(i,j) = interp1(t,ksi_y_dot(:,j),interp_t(i));
    end
    j
end

% Вычисление вторых производных
dt = mean(diff(interp_t)); 
for j=1:length(x_dot_inrp(1,:))
   for i=1:length(x_dot_inrp(:,1))
       if i == 1
           x_dd(i,j) = (x_dot_inrp(i+1,j)-x_dot_inrp(i+1,j))/dt;
           y_dd(i,j) = (y_dot_inrp(i+1,j)-y_dot_inrp(i+1,j))/dt;
       elseif i == length(x_dot_inrp(:,1))
           x_dd(i,j) = (x_dot_inrp(i-1,j)-x_dot_inrp(i,j))/dt;
           y_dd(i,j) = (y_dot_inrp(i-1,j)-y_dot_inrp(i,j))/dt;
       else
           x_dd(i,j) = (x_dot_inrp(i+1,j)-x_dot_inrp(i-1,j))/(2*dt);
           y_dd(i,j) = (y_dot_inrp(i+1,j)-y_dot_inrp(i-1,j))/(2*dt);
       end 
   end
end

% Вычисления момента в каждом сечении
h = 1/m;
epsilon = @(x) (1e-3/0.7)*sin(pi*x);
pfi = @(x) 1*sin(2*pi*x);
absciss = 0:h:1;
for i = 1:length(absciss)
    eps_v(i) = epsilon(absciss(i));
    pfi_v(i) = pfi(absciss(i));
end

for j=1:length(pfi_v)       % столбцы (сечения)
   for i=1:length(interp_t) % строки (время)
       det1 = (1/m)*(x_inrp(i,j)*(y_dd(i,j)-N^2*eps_v(j)*sin(pfi(j)+N*interp_t(i)))+eps_v(j)*cos(pfi(j)+N*interp_t(i))*(y_dd(i,j)-N^2*eps_v(j)*sin(pfi(j)+N*interp_t(i))));
       det2 = (1/m)*(x_dd(i,j)*(y_inrp(i,j)+eps_v(j)*sin(pfi(j)+N*interp_t(i)))-N^2*eps_v(j)*cos(pfi(j)+N*interp_t(i))*(y_inrp(i,j)+eps_v(j)*sin(pfi(j)+N*interp_t(i))));
       Moment(i,j) = det1-det2;
   end
end

% Вычиследние момента в каждый момент времени
for i=1:length(Moment(:,1))
    Sum_m(i) = sum(Moment(i,:));
end


end