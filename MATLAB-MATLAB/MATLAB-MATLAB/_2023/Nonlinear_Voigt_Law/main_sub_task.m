close all
clc
clear
% Кубический закон Фойхта
% Схема: заделка-заделка
warning('on')
m = 6;
h = 1/m;
N = 2;
Mom = 2;

% Параметры стержня
d = 20*10^-3;           %[m]
l = 0.7;                %[m]
ro = 7800;              %[kg/m^3]
Elastic = 2e11;         %[Pa]
zeta_e = 0.025;         %[-]
zeta_V = 0.005;         %[-]

J = pi*d^4/64;          %[m^4]
C = Elastic*J;
eps_d = d/l;            %[-]
beta_R = (eps_d^2)/16;  %[-]
mu_R = 1;

%% Получение матриц
[G00Int,G01Int,G02Int,G20Int,G21Int,...
          G00,G01,G02,G20,G21,J]=get_matrix(m);

%%
I=eye(m+1);         % единичная матрица
E = eye(2);         % единичная матрица
EC = zeros(m+1); EC(end,end) = 1;
Z=zeros(m+1);       % нулевая матрица
R = [0 -1; 1 0];    

%% Исключаем первый и последний узлы
[G00Int,G01Int,G02Int,G20Int,G21Int,...
          G00,G01,G02,G20,G21,J,I,EC,Z]=matrix_edging(G00Int,G01Int,G02Int,G20Int,G21Int,...
          G00,G01,G02,G20,G21,J,I,EC,Z,m);
      
%% Заполнение матриц коэффициентов

a_ksi_0 = kron(I,(E-2*zeta_V*N*R));
a_theta_0 = Mom*kron(G01Int,E);

b_ksi_0 = kron(Z,E);
b_theta_0 = kron(I,E-2*zeta_V*N*R)+Mom*kron(G21Int,R);

A0 = [a_ksi_0, a_theta_0;
      b_ksi_0, b_theta_0];

a_ksi_1 = 2*zeta_V*kron(I,E)+2*zeta_e*kron(G00Int,E);
a_theta_1 = 2*beta_R*N*kron(G00Int,E);

b_ksi_1 = 2*zeta_e*kron(G20Int,R);
b_theta_1 = 2*zeta_V*kron(I,E)+2*beta_R*N*kron(G20Int,R);

A1 = [a_ksi_1, a_theta_1;
      b_ksi_1, b_theta_1];

  
a_ksi_2 = mu_R*kron(G00Int,E);
a_theta_2 = beta_R*kron(G00Int,R);

b_ksi_2 = mu_R*kron(G20Int,R);
b_theta_2 = -beta_R*kron(G20Int,E);


A2 = [a_ksi_2, a_theta_2;
      b_ksi_2, b_theta_2];




% g = [kron(G00Int,E), zeros(length(kron(G00Int,E)));zeros(2,length(kron(G00Int,E))) ,kron(G30Int(end,:),E)];


%% Поиск собственных значений
EigVal = polyeig(A0,A1,A2);

figure;
hold on; grid on; box on;
ff = gca;
ff.FontSize = 20;
ff.FontName = 'Times New Roman';

flagNonStab = 0;
for i=1:length(EigVal)
    markerColor = '.k';
    markerSize = 16;
    if real(EigVal(i))>0
        markerColor = '.r';
        markerSize = 25;
        flagNonStab = flagNonStab+1;
    end
    plot(real(EigVal(i)),imag(EigVal(i)),markerColor,'MarkerSize',markerSize)
end
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')

if flagNonStab>0
    title(['Система НЕустойчива ','(max(Re \lambda) = ',num2str(max(real(EigVal))),')'])
else
    title(['Система усточива ','(max(Re \lambda) = ',num2str(max(real(EigVal))),')'])
end

%% Поиск критической скорости
N_start = 0;
step_N = 0.5;

Argan_diagram (N_start,step_N);


%% Решение ДУ

opt=odeset('AbsTol',1e-7,'RelTol',1e-7); % настройки точности

T=[0,1000];   % интервал полного исследования
x0 = zeros(2*length(A0),1); % Вектор начальных условий
amplitude = 0;
period = 2;
phase = 1;
tic;
[t,Z0] = ode15s(@(t,Z0) solver(t,Z0,A0,A1,A2,g,h,N,kappa_B, eta_B,amplitude,phase,period),T,x0,opt);
toc;


%% Формирование ksi_x и ksi_y в каждом узле
flag = 1;
for j=1:2:(length(A2)-2)
    for i=1:length(Z0(:,1))
        ksi_x(i,j) = Z0(i,j) + Z0(i,length(A0)-1)*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_x(:,2:2:end) = [];
ksi_x = [zeros(length(ksi_x(:,1)),1),ksi_x,Z0(:,length(A0)-1)];

flag = 1;
for j=2:2:(length(A2)-2)
    for i=1:length(Z0(:,1))
        ksi_y(i,j) = Z0(i,j) + Z0(i,length(A0))*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_y(:,1:2:end) = [];
ksi_y = [zeros(length(ksi_y(:,1)),1),ksi_y,Z0(:,length(A0))];

z_vect = 0:h:1;

%% Формирование скоростей ksi_x и ksi_y в каждом узле
flag = 1;
for j=length(A2)+1:2:2*(length(A2))-2
    for i=1:length(Z0(:,1))
        ksi_x_dot(i,j) = Z0(i,j) + Z0(i,2*length(A0)-1)*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_x_dot(:,1:length(A2))= [];
ksi_x_dot(:,2:2:end) = [];
ksi_x_dot = [zeros(length(ksi_x_dot(:,1)),1),ksi_x_dot,Z0(:,2*length(A0)-1)];

flag = 1;
for j=length(A2)+2:2:2*(length(A2))-1
    for i=1:length(Z0(:,1))
        ksi_y_dot(i,j) = Z0(i,j) + Z0(i,2*length(A0))*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_y_dot(:,1:length(A2)+1)= [];
ksi_y_dot(:,2:2:end) = [];
ksi_y_dot = [zeros(length(ksi_y_dot(:,1)),1),ksi_y_dot,Z0(:,2*length(A0))];

z_vect = 0:h:1;

%%
figure;
box on; grid on; hold on;
plot(t,ksi_x(:,end))


%% 
t_t_vect = linspace(t(end)/1.5,t(end),10000);
ksi_0_x = @(t_t) (amplitude*sin((2*pi/period)*t_t+phase));
ksi_0_y = @(t_t) (amplitude*cos((2*pi/period)*t_t+phase));

for i=1:length(t_t_vect)
    ksi_0_x_(i) = ksi_0_x(t_t_vect(i));
    ksi_0_y_(i) = ksi_0_y(t_t_vect(i));
    ksi_X_1(i) = interp1(t,ksi_x(:,end),t_t_vect(i));
    ksi_Y_1(i) = interp1(t,ksi_y(:,end),t_t_vect(i));
end
%%
figure1 = figure;
axes1 = axes('Parent',figure1);
set(axes1,'DataAspectRatio',[1 1 1]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'on');

plot(ksi_0_x_,ksi_0_y_)
plot(ksi_X_1,ksi_Y_1)

xlabel('\xi_{x}')
ylabel('\xi_{y}','Rotation',0)
legend('Support','Beam')
ff = gca;
ff.FontSize = 16;

%% Трехмерный график перемещений
for i=1:length(z_vect)
    Z_3D{i} = kron(ones(length(ksi_x),1),z_vect(i));
end
index_for_time = find(t<=t_t_vect(end)&t>=t_t_vect(1));
figure;
hold on; box on; grid on;
view([45.9000000862685 13.7999998741406])
for i=1:length(z_vect)
    plot3(Z_3D{i}(index_for_time),ksi_x(index_for_time,i),ksi_y(index_for_time,i),'.-k')
end

%% Трехмерный график 2
% figure;
% hold on; box on; grid on;
% view([34.8750002358158 46.0463697962445])
% for i=1:length(index_for_time)
%     plot3(z_vect,ksi_x(i,:),ksi_y(i,:),'.-b')
%     drawnow 
%     pause(0.1)
% end
% 
% ff = gca;
% ff. FontSize = 16; 
% xlabel('z')
% ylabel('x')
% zlabel('y')
pause(1);
%%
[t_m,Mi] = get_moment(t,ksi_x,ksi_y,ksi_x_dot,ksi_y_dot,m,N);
%%
t_m(end) = [];
Mi(end) = [];
%%
figure;
box on; grid on; hold on
plot(t_m,Mi);
xlabel('time')
ylabel('M')
ff = gca;
ff.FontSize = 16;

