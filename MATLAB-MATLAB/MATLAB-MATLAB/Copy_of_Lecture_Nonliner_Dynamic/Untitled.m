close all
clc
clear
%%
% Глобальные переменные
global I E S Z eta_i eta_e eta_Dksi eta_Dte
global G00 G10 G01 G11 G00Int G01Int G10Int G11Int
global DE EC beta muR betaR h

%%%%%%%%%%%%%%%%%%%%%%%%__ЗАДАНИЕ ПАРАМЕТРОВ__%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Размерные параметры
ro = 7800;                  % [кг/м^3] - плотность материала
Elastic = 2.1*10^11;        % [Па] - модуль упругости материала 
l = 0.7;                    % [м] - длина стержня
d = 20*10^-3;               % [м] - диаметр стержня 
D = 0.2;                    % [м] - диаметр диска
a = 10*10^-3;               % [м] - толщина диска

h_j = 10*10^-3;             % [м] - зазор между диском и ограничителями
e = 1*10^-3;                % [м] - величина эксцентриситета
gravity = 0;             % [м/c^2] - ускорение свободного падения

% Безразмерные параметры
eta_i = 0.005;              % [-] - коэффициент внутреннего линейного демпфирования стержня
eta_e = 0.025;              % [-] - коэффициент внешнего линейного демпфирования стержня
eta_Dksi = 0.0025;          % [-] - коэффициент внешнего линейного демфирования диска 
eta_Dte = 0.05;             % [-] - коэффициент внешнего углового демфирования диска

kappa_j = 850;              % [-] - коэффициент жесткости ограничителей
d_j = 0.06;                 % [-] - коэффициент демпфирования ограничителей   
f_j = 0.1;                  % [-] - коэффициент сухого трения c ограничителями

N = 16;                     % [-] - cкорость вращения ротора
m = 6;                     % количество участков разбиения
%NumbOp=60:120:300;          %[degree] - расположение опор (от Ox против часовой стрелки)
NumbOp=180;
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
%% Подготавливаем матрицы нужных размерностей
L=m+1;        % количество узлов
E=eye(2);     % единачная матрица 2x2
S=[0 1;-1 0]; % кососимметричная матрица 2x2
I=eye(L);     % единичная матрица (кол-во узлов)x(кол-во узлов)
Z=zeros(L);   % нулевая матрица (кол-во узлов)x(кол-во узлов)

EC=zeros(L);       % нулевая матрица (кол-во узлов)x(кол-во узлов)
EC(jC,jC) = 1;     % единственный ненулевой – последний диагональный элемент

E1=zeros(L);       % нулевая матрица (кол-во узлов)x(кол-во узлов)
E1(1,1)=1;         % единственный ненулевой – 1 диагональный элемент

EL=zeros(L);       % нулевая матрица (кол-во узлов)x(кол-во узлов)
EL(L,L)=1;         % единственный ненулевой – последний диагональный элемент

DE=EL-E1;          % матрица разности 
%% --------------------------
% Подготавливаем матрицы Грина (пока все нулевые)
G00=zeros(L);G01=zeros(L);
G10=zeros(L);G11=zeros(L);
G00Int=zeros(L);G01Int=zeros(L);
G10Int=zeros(L);G11Int=zeros(L);

% Заполняем матрицы Грина

for j=1:L     % проходим по строкам матриц Грина
    for k=1:L % проходим по столбцам матриц Грина
     [G00(j,k),G01(j,k),G10(j,k),G11(j,k)]=MatrixOfGreen(j,k,h,System_name);
     [G00Int(j,k),G01Int(j,k),G10Int(j,k),G11Int(j,k)]=MatrixOfGreenInt(j,k,h,System_name);
    end
end
%% Окаймление матриц – убираем заведомо однородные степени свободы
[G00,G01,G10,G11,G00Int,G01Int,G10Int,G11Int,I,Z,EC,m,jC,DE,mnew] = Matrix_Edging(G00,G01,G10,G11,G00Int,G01Int,G10Int,G11Int,I,Z,EC,DE,m,jC,L,System_name);

%% Заполнение матриц коэффициентов
aksi=kron(I,E+2*eta_i*N*S);
dksi=kron(Z,E);
ateta=kron(Z,E);
dteta=kron(I,S-2*eta_i*N*E);

bksi=kron(2*eta_i*I+2*eta_e*h*G00Int+2*eta_Dksi*G00*EC,E);
eksi=kron(2*betaR*N*(-h*G01Int+G00*DE)-2*beta*N*G01*EC,E)+kron(2*eta_Dte*G01*EC,S); 
bteta=kron(2*eta_e*h*G10Int+2*eta_Dksi*G10*EC,E);
eteta=kron(2*betaR*N*(-h*G11Int+G10*DE)-2*beta*N*G11*EC,E)+kron(2*eta_i*I+2*eta_Dte*G11*EC,S);


cksi=kron(muR*h*G00Int+G00*EC,E);
fksi=kron(betaR*(h*G01Int-G00*DE)+beta*G01*EC,S);
cteta=kron(muR*h*G10Int+G10*EC,E);
fteta=kron(betaR*(h*G11Int-G10*DE)+beta*G11*EC,S);

A0=[aksi dksi;ateta dteta];

A1=[bksi eksi;bteta eteta];

A2=[cksi fksi; cteta fteta];

% Матрица g коэффициентов правой части
ZeroMatr = kron(Z,E);
gksi = kron(G00*EC,E);
gteta = kron(G10*EC,E);
g = [gksi,ZeroMatr;ZeroMatr,gteta];
% Матрица h коэффициентов правой части
hksi = kron(G00Int,E);
hteta = kron(G10Int,E);
h_matr = [hksi,ZeroMatr;ZeroMatr,hteta];
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
    tNS = 30; % Количство последних оборотов
    ind = 35;
Moving_the_center(time{ind},ZZ{ind},TN,tNS,disp,angle,N)

%% Траектория центра диска со стробоскопическим отображением
    tNS = 200; % Количство последних оборотов
    ind = 35;




tN=time{ind}/TN;                                                % переходим от безразмерного времени к количеству оборотов
tPoincare=(time{ind}(1):TN:time{ind}(end))';                            % вектор времени, разбитый на равные интервалы, равные времени одного оборота ротора
ZxPoincare=interp1(time{ind},ZZ{ind}(:,disp(1)),tPoincare,'spline'); % вектор перемещений по oX через каждый оборот ротора 
ZyPoincare=interp1(time{ind},ZZ{ind}(:,disp(2)),tPoincare,'spline'); % вектор перемещений по oY через каждый оборот ротора
tPoincareIndex=tPoincare/TN;                            % вектор количества оборотов       
II=find(tPoincareIndex>(tPoincareIndex(end)-tNS));     % достаём индексы точек стробирования, соотвестствующие tNS2 последним оборотам
I=find(tN>(tN(end)-tNS));                               % достаём индексы точек самой траекториии, соотв. tNS последним оборотам                                


% Построение фазовой траектории с точкам стробирования
X1 = ZZ{ind}(I,disp(2));
Y1 = ZZ{ind}(I,disp(1));
X1_poin = ZyPoincare(II);
Y1_poin = ZxPoincare(II);
figure('WindowState','maximized'); 
    hold on
    plot(-X1,Y1,'LineWidth',0.5);
    plot(-X1_poin,Y1_poin,'r.','MarkerSize',18);
    axis equal
    xlabel('\xi_{\it y}')
    ylabel('\xi_{\it x}','Rotation',0)
    %title(['N = ',num2str(N),'; Количество последних оборотов: ',num2str(tNS)])
    ax1 = gca;
    ax1.FontName = 'Times New Roman';
    ax1.FontSize = 18;
    grid on; hold on; box on;
    xlim padded
    ylim padded
    
%%
figure('WindowState','maximized'); 
    hold on
    plot(-X1,Y1,'LineWidth',0.5);
    plot(-X1_poin,Y1_poin,'r.','MarkerSize',18);
    
    plot(-X2,Y2,'LineWidth',0.5);
    plot(-X2_poin,Y2_poin,'r.','MarkerSize',18);
    
    plot(-X3,Y3,'LineWidth',0.5);
    plot(-X3_poin,Y3_poin,'r.','MarkerSize',18);
    axis equal
    xlabel('\xi_{\it y}')
    ylabel('\xi_{\it x}','Rotation',0)
    %title(['N = ',num2str(N),'; Количество последних оборотов: ',num2str(tNS)])
    ax1 = gca;
    ax1.FontName = 'Times New Roman';
    ax1.FontSize = 18;
    grid on; hold on; box on;
    xlim padded
    ylim padded

%%
figure('WindowState','maximized'); 
    hold on
    plot(-X3,Y3,'LineWidth',0.5,'Color','#7E2F8E');
    plot(-X3_poin,Y3_poin,'r.','MarkerSize',18);
    
    axis equal
    xlabel('\xi_{\it y}')
    ylabel('\xi_{\it x}','Rotation',0)
    %title(['N = ',num2str(N),'; Количество последних оборотов: ',num2str(tNS)])
    ax1 = gca;
    ax1.FontName = 'Times New Roman';
    ax1.FontSize = 18;
    grid on; hold on; box on;
    xlim padded
    ylim padded







