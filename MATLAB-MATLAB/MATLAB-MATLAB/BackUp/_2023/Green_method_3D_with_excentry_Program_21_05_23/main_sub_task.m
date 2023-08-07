close all
clc
clear
warning('on')

%% Задание параметров %%

% Конфигурация задачи
m = 30;                 % Кол-во участков
N = 2;                  %[1] Скорость вращения
Mom = 2;                %[1] Крутящий момент

% Параметры стержня
d = 20*10^-3;           %[m]
l = 0.7;                %[m]
ro = 7800;              %[kg/m^3]
Elastic = 2e11;         %[Pa]
zeta_e = 0.025;         %[-]
zeta_V = 0.005;         %[-]
Ampl_eps = 1e-3;        %[m]    Амплитуда дисбаланса
Ampl_phase = 1;         %[rad]  Амплитуда фазы



% Параметры концевой массы
m_B = 0.0001;                %[кг]
k_B = {[],0};         %[N/m](задайте размерную величину или долю от изг.ж)
eta_B = 0.1;            %[-]

% Параметры кинематического нагружения
amplitude = 0;          %[1]
period = 2;             %[1]
phase = 1;              %[rad]

% Параметры решения
opt=odeset('AbsTol',1e-7,'RelTol',1e-7);
nT = 500;              % количество оборотов ротора
initial_vect = [];     % если пустой, то нулевой


%%%%%%%%%%%%%%%%%%%%%%%%

%% Обработка данных
[dat,sol] = get_struct(m,N,Mom,d,l,ro,Elastic,zeta_e,zeta_V,m_B,k_B,eta_B,...
                   Ampl_eps, Ampl_phase,...
                   amplitude, period, phase,nT,initial_vect);


%% Получение матриц
[dat] = get_matrix(dat);

%% Исключение первого узла
[dat] = fringing(dat);

%% Получение матриц-коэффициентов
[dat] = get_matr_coeff(dat);

%% Поиск собственных значений при выбиранном N
get_eigenValues(dat);

%% Поиск критической скорости (ваирование скорости)
N_start = 0;    % Начальная точка
step_N = 0.1;  % Шаг поиска
try
    Argan_diagram (N_start,step_N,dat);
    
catch exception
    disp('Критическую скорость не удалось найти :(')
end


%% Решение ДУ
tic;
[t,Z0] = ode23t(@(t,Z0) solver(t,Z0,dat),sol.T,sol.x0,opt);
toc;

%% Формирование ksi_x, ksi_y и их скоростей
[ksi_x, ksi_y,ksi_x_dot,ksi_y_dot] = get_KSI(dat,Z0); 

%% Вычисление момента
% [t_m,Mi] = get_moment(t,ksi_x,ksi_y,ksi_x_dot,ksi_y_dot,dat);
%%
                %%% Визуализация %%%

                              
%% Графики перемещений за N_end - последних оборотов
N_end = 250;

disp_plot(N_end,t,ksi_x,ksi_y,dat,sol)

%% Трехмерный график центров сечений
%       за N_end - последних оборотов
N_end = 250;

plot_3D(N_end,t,ksi_x,ksi_y,dat,sol)

%% График крутящего момента за N_end последних оборотов
N_end = 250;
% moment_plot(N_end,t_m,Mi,dat,sol)

%% Информация о расчёте
information_plot(dat)
%%
titlE=['\beta = ',          num2str(beta),...
          ', \zeta_{ i} = ',   num2str(zeta_i),...
          ', \zeta_{ e} = ',   num2str(zeta_e),...
          ', Z_{ e} = ',       num2str(Zeta_e),...
          ', Z_{ \theta e} = ',num2str(Zeta_te),...
          ', \zeta_{ C} = ',   num2str(z_C)];

