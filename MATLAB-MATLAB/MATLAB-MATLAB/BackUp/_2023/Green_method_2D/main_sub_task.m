close all
clc
clear

m = 7;
h = 1/m;
N = 22;

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
%beta = (eps_d^2)/16;    %[-]
beta = 0;
% Параметры краевого условия
m_B = 1;               %[кг]
k_B = 400;              %[N/m];
%d_B = ///              %

mu_B = m_B/(ro*(pi*d^2/4)*l*l);
kappa_B = k_B*l^3/C;
kappa_B = 800;
%eta_B = d_B*l/(2*sqrt(ro*(pi*d^2/4)*l*C));
eta_B = 0.001;          %[-]
eta_B = 1;

%% Получение матриц

[Const_F1,Const_F2,Const_F1_ticks,Const_F2_ticks,T,G00Int,G02Int,G30Int]=get_matrix(m);

I=eye(m+1);         % единичная матрица
Z=zeros(m+1);       % нулевая матрица
Z1=zeros(m+1,1);    % нулевой вектор

%% Исключаем первый и послений узлы
Const_F1 = Const_F1(2:m);
Const_F2 = Const_F2(2:m);
T = T(2:m);
G00Int = G00Int(2:m,2:m);
G02Int = G02Int(2:m,2:m);
G30Int = G30Int(2:m,2:m);
I = I(2:m,2:m);
Z = Z(2:m);
Z1 = Z1(2:m);


%% Основные матрицы

A0 = [I,    Z1;...
      Z1',  12+kappa_B];

A1 = [2*zeta_V*I+2*zeta_e*G00Int,   2*zeta_e*Const_F2;...
      2*zeta_e*G30Int(end,:),       2*zeta_e*Const_F2_ticks+24*zeta_V+2*eta_B];

A2 = [G00Int-beta*G02Int            , -beta*Const_F1+Const_F2;...
      G30Int(end,:)-beta*T   , mu_B-beta*Const_F1_ticks+Const_F2_ticks];

%%



opt=odeset('AbsTol',1e-8,'RelTol',1e-8); % настройки точности

T=[0,100];   % интервал полного исследования
x0 = zeros(2*length(A0),1); % Вектор начальных условий
amplitude = 0.05;
period = 0.7007;
phase = 1;
tic;
x0 = zeros(2*length(A0),1);

[t,Z0] = ode23t(@(t,Z0) solver(t,Z0,A0,A1,A2,kappa_B, eta_B,amplitude,phase,period),T,x0,opt);
toc;

                      %%% Визаулизация %%%
%% Формирование ksi_x в каждом узле
ksi_x = zeros(length(Z0(:,1)),(length(A2)-1));
for j=1:(length(A2)-1)
    for i=1:length(Z0(:,1))
        ksi_x(i,j) = Z0(i,j) + Z0(i,length(A0))*(3*(j*h)^2-2*(j*h)^3);
        z_vect{j} = j*h*ones(length(Z0(:,1)),1);
    end
end
%%
z_vect{end+1} = ones(length(Z0(:,1)),1); 
ksi_x = [ksi_x,Z0(:,(length(A2)))];
% %% Для всех времен
% figure;
% hold on; box on; grid on;
% for i=1:length(ksi_x(1,:))
%     plot(z_vect{i}(:),ksi_x(:,i));
% end 
% ff = gca;
% ff.FontSize = 15;
% xlabel('координата \zeta')
% ylabel('\xi_{ \itx}')
% %% Векторы для построения
% for ttime=1:length(ksi_x(:,1))
%     KSI{ttime} = ksi_x(ttime,:);
%     KSI{ttime} = [0,KSI{ttime}];
% end
 for ttime=1:length(ksi_x(1,:))
     ZETA(ttime) = z_vect{ttime}(1);
end
ZETA=[0,ZETA];
% %% Для определенного момента времени
% figure;
% hold on; box on; grid on;
% plot(ZETA,KSI{492},'.-r','MarkerSize',15);

%% Интерполяционные значения
t_interp = 90;
for i=1:length(ksi_x(1,:))
    ksi_x_interp(i) = interp1(t,ksi_x(:,i),t_interp);
end
ksi_x_interp = [0,ksi_x_interp];
%%
figure;
hold on; box on; grid on;
plot(ZETA,ksi_x_interp,'.-r','MarkerSize',15);
ff = gca; 
ff.FontSize = 18;
xlabel('\zeta');
ylabel('\xi_{ x}');
txt = ['Для t = ', num2str(t_interp),'; phase = ', num2str(phase),'; ampl = ' num2str(amplitude)];
title(txt)


%% Сравнение с нагружением

Ksi_x_point = @(x) interp1(t,ksi_x(:,end),x);

Ksi_0 = @(x) amplitude*sin((2*pi/period)*x+phase);

x_vect = linspace(0,t(end),1e8);
%%
figure;
grid on; hold on; box on;
plot(x_vect,Ksi_x_point(x_vect),'-r','LineWidth',1.5)
plot(x_vect,Ksi_0(x_vect),'-b','LineWidth',1.5) 
legend('\xi_{ x}','\xi_{ 0}')
ff = gca;
ff.FontSize = 20;


