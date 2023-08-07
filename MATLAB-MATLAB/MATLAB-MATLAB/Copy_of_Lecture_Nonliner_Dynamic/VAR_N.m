close all
clc
clear
%%
m = 6; 
mnew = m-1;
x0 = zeros(8*(mnew+1),1);
% N_vect = 11:0.25:17;
N_vect = 11:0.25:13;
NUMREV = 300;
%%
for i = 1:length(N_vect)
    disp(i);
   [t,Z0,XextrEta,YextrEta,x_initial,Average]=for_var_N(N_vect(i),x0,NUMREV);
    x0 = x_initial;
    time{i} = t;
    Z{i} = Z0;
    EXTRx{i} = XextrEta; EXTRy{i} = YextrEta;
    AVGL{i} = Average;
end

%% Построение бифуркационной диаграммы
%%
for i = 1:1:length(AVGL)
   avgL(i) = AVGL{i};
end
N_crit = 13.9185;
figure;
subplot(211)
hold on; box on; grid on;
for k = 1:length(N_vect)
    marker = '.k';
    if AVGL{k}<0
        marker = '.r';
    end 
    plot(N_vect(k),EXTRx{k},marker,'MarkerSize',18);
end
    ff = gca; 
    ff.FontName = 'Times New Roman';
    ff.FontSize = 20;
    xlabel('N');
    ylabel('EXTR[\xi_{\it x}]')
    xlim([10.75,17.25])
    ylim([-0.03,0.03])
subplot(212)
hold on; box on; grid on;
% for k = 1:length(N_vect)
%     plot(N_vect(k),avgL(k),'.r','MarkerSize',18)
% end
plot(N_vect,avgL,'.-r','MarkerSize',22,'LineWidth',2,'Color','#7E2F8E')
xline(N_crit,'--r',[{'N_{crit}'};'(',num2str(N_crit),')'],'LabelVerticalAlignment','bottom','LineWidth',2,'FontName','Times New Roman','FontSize',14);
yline(0,'--k')
    ff = gca;
    ff.FontName = 'Times New Roman';
    ff.FontSize = 20;
    xlabel('N');
    ylabel('Avg(\Lambda)')
    xlim([10.75,17.25])
    ylim([-1,1.5])

% %% Совместные графики
% figure; 
% hold on; box on; grid on;
% plot(N_vect,avgL1,'.-','MarkerSize',22,'LineWidth',2,'Color','#0072BD')
% plot(N_vect,avgL2,'.-','MarkerSize',22,'LineWidth',2,'Color','#D95319')
% plot(N_vect,avgL3,'.-','MarkerSize',22,'LineWidth',2,'Color','#7E2F8E')
% yline(0,'--k','LineWidth',2)
% xline(N_crit,'--r',[{'N_{crit}'};'(',num2str(N_crit),')'],'LabelVerticalAlignment','bottom','LineWidth',2,'FontName','Times New Roman','FontSize',18);
%     ff = gca;
%     ff.FontName = 'Times New Roman';
%     ff.FontSize = 20;
%     xlabel('N');
%     ylabel('Avg(\Lambda)')   
%     ylim([-1,1.5]);
%     xlim('padded')
%     legend('Одна опора','Две опоры','Три опоры')
    
%%
% Размерные параметры
ro = 7800;                  % [кг/м^3] - плотность материала
Elastic = 2.1*10^11;        % [Па] - модуль упругости материала 
l = 0.7;                    % [м] - длина стержня
d = 20*10^-3;               % [м] - диаметр стержня 
D = 0.2;                    % [м] - диаметр диска
a = 10*10^-3;               % [м] - толщина диска

h_j = 10*10^-3;             % [м] - зазор между диском и ограничителями
e = 1*10^-3;                % [м] - величина эксцентриситета
gravity = 9.81;             % [м/c^2] - ускорение свободного падения

% Безразмерные параметры
eta_i = 0.005;              % [-] - коэффициент внутреннего линейного демпфирования стержня
eta_e = 0.025;              % [-] - коэффициент внешнего линейного демпфирования стержня
eta_Dksi = 0.0025;          % [-] - коэффициент внешнего линейного демфирования диска 
eta_Dte = 0.05;             % [-] - коэффициент внешнего углового демфирования диска

kappa_j = 850;              % [-] - коэффициент жесткости ограничителей
d_j = 0.06;                 % [-] - коэффициент демпфирования ограничителей   
f_j = 0.1;                  % [-] - коэффициент сухого трения c ограничителями
                            % [-] - cкорость вращения ротора
m = 6;                     % количество участков разбиения
NumbOp=60:120:300;          %[degree] - расположение опор (от Ox против часовой стрелки)
zDisk=1;                  %[-] - кооридната расположения диска (*в пределах (0;1)*)
System_name=2;              % (1) "заделка – заделка"; (2) "консоль"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Расчёт параметров
M = ro*(pi*D^2/4)*a;                % [кг] - масса диска
B = M*D^2/16;                       % [кг*м^2] - физический момент инерции диска
IR = pi*d^4/64;                     % [м^4] - геометрический момент инерции сечения стержня
mR = ro*pi*d^2/4;                   % [кг/м] - погонная масса стержня
beta = B/(M*l^2);                   % безразмерный коэффициент момента инерции диска
betaR = ro*IR/(M*l);                % безразмерный коэффициент поворотной инерции 
muR = mR*l/M;                       % безразмерный коэффициент массы стержня
chi_j = h_j/l;                      % коэффициент зазора
eD = D/l;                           % коэффициент диаметра диска
epsilon = e/l;                      % коэффициент эксцентриситета
gamma = gravity*M*l^2/(Elastic*IR); % коэффициент силы тяжести

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=1/m; % шаг
% Граничные случаи задание постановки диска
s=round(zDisk/h);
%if s==0         % Если поставить диск в левую заделку
  %  jC=2;       
 %elseif s==m    % Если поставить диск в правую заделку
 %    jC=m;    
%else
    jC=s+1; % нужно смещение вправо на один узел
%end
if ((zDisk>=1||zDisk<=0)&&System_name==1)
    error ('Проверьте постановку диска');
end
if ((zDisk~=1)&&System_name==2)
    error ('Проверьте постановку диска');
end
jC = jC - 1;
%% Описание индексации задачи

index_disp = 1:(mnew+1)*2;                  % Степени свободы перемещений

disp = [index_disp(jC*2)-1;...              % ksi_x где расположен диск
        index_disp(2*jC)];                  % ksi_y где расположен диск 

index_angle = (mnew+1)*2+1:(mnew+1)*4;      % Степени свободы углов поворота

angle = [index_angle(jC*2)-1;...            % teta_x где расположен диск
         index_angle(2*jC)];                % teta_y где расположен диск 

% То же самое, только для скоростей     

index_disp_vel = (mnew+1)*4+1:(mnew+1)*6;   

disp_vel = [index_disp_vel(jC*2)-1;...  
            index_disp_vel(2*jC)];
   
index_angle_vel = (mnew+1)*6+1:(mnew+1)*8;

angle_vel = [index_angle_vel(jC*2)-1;index_angle_vel(2*jC)];

%% Графики перемещений и углов поворота
    tNS = 10; % Количство последних оборотов
    ind = 9;
    TN=2*pi/N_vect(ind);
    
    Moving_the_center(time{ind},Z{ind},TN,tNS,disp,angle,N_vect(ind))

%% Траектория центра диска со стробоскопическим отображением
    tNS = 200; % Количство последних оборотов
    ind = 9;
Phase_trajectory_2D(time{ind},Z{ind},TN,tNS,disp,N_vect(ind))

%% Квазифазовая траектория центра диска
    tNS = 300; % Количство последних оборотов
    ind = 9;
Quasi_phase_trajectorys(time{ind},Z{ind},TN,tNS,disp)

%% Спектральный анализ перемещений
    tNS = 200; % Количество последних оборотов
    ind = 9;
Spectral_analysis(time{ind},Z{ind},tNS,TN,disp)

%% Реакции в упругодемпфирующих опорах
    tNS = 20; % Количство последних оборотов
    ind = 9;
Reactions(time{ind},Z{ind},TN,tNS,NumbOp,disp,disp_vel,kappa_j,h_j,d_j,f_j,N_vect(ind),eD)
%% Коэффициент прецессии
    tNS = 200; % Количство последних оборотов
    ind = 9;
Precession_coefficient(time{ind},Z{ind},TN,tNS,disp,N_vect(ind))

%% Построение положение узлов в 3D
    tNS = 50; % Количство последних оборотов
    ind = 9;
plot_Disp_3D(time{ind},Z{ind},TN,tNS,m,index_disp,jC,System_name)


