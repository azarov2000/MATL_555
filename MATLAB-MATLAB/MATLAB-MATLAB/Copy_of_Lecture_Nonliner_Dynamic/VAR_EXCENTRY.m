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
ro = 7800;                      % [кг/м^3] - плотность материала
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

N = 16;                     % [-] - cкорость вращения ротора
m = 6;                     % количество участков разбиения
%NumbOp=60:120:300;          %[degree] - расположение опор (от Ox против часовой стрелки)
NumbOp = 0:180:180;
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

%% Построение диаграммы Аргана
N_initial = 0; % начальная точка
step_N = 0.1; % Шаг изменения скорости
        
Argan_diagram (N_initial,step_N)

%%
%%%%%%%%%%%%%%%%%%%%%%____РЕШЕНИЕ СИСТЕМЫ____%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt=odeset('AbsTol',1e-6,'RelTol',1e-6); % настройки точности
TN=2*pi/N;    % время одного оборота
nT = 6000;    % количество оборотов ротора
tlim=TN*nT;   % полное время исследования
T=[0,tlim];   % интервал полного исследования
x0 = zeros(8*(mnew+1),1);  x0(disp(:)) = 0.001*[1,1]; % Вектор начальных условий

%% **************************___РЕЗУЛЬТАТЫ___******************************
%% Построение корней хар. уравнения для при всех выбранных параметрах

% %% ДОПОЛНЕНИЕ (получение биффуркационной диаграммы и запоминание векторов состояния при различных параметрах)
vect_ex = linspace(0,2*epsilon,40);
le_vect_ex = length(vect_ex);
NUMREV=300; % количество последних оборотов
for j=1:1:le_vect_ex
    j 
    [t1,Z1]=ode23t(@(t,Z0) rGap_Multiple_Supports(t,Z0,A0,A1,A2,g,h_matr,N,eD,mnew,muR,chi_j,kappa_j,d_j,vect_ex(j),gamma,f_j,NumbOp,disp,disp_vel),T,x0,opt);
    tEta=t1/TN;
    NumberOfRev=find(tEta>(tEta(end)-NUMREV));
    [textrEtaX,XextrEta] = ext(tEta(NumberOfRev),Z1(NumberOfRev,disp(1)));
    [textrEtaY,YextrEta] = ext(tEta(NumberOfRev),Z1(NumberOfRev,disp(2)));
    MaxKSIx{j}=XextrEta;
    MaxKSIy{j}=YextrEta;
    % Сохранение реализаций
    ZZ{j} = Z1; time{j} = t1;
    x0 = Z1(end,:);
    
    tay{j} = tEta(NumberOfRev)*TN;
    
    Px=polyfit(tay{j},Z1(NumberOfRev,disp(1)),0);
    Py=polyfit(tay{j},Z1(NumberOfRev,disp(2)),0);
    TrandX=polyval(Px,tay{j});    % Получаем записываем ординаты этого полинома
    TrandY=polyval(Py,tay{j});    % Получаем записываем ординаты этого полинома
 
    Z1(NumberOfRev,disp(1))=Z1(NumberOfRev,disp(1))-TrandX; % по сути мы убираем постоянную составляющую сигнала
    Z1(NumberOfRev,disp(2))=Z1(NumberOfRev,disp(2))-TrandY;
    
    [tPr,CoeffPr] = PrecessionCoeff(tay{j},Z1(NumberOfRev,disp(1)),Z1(NumberOfRev,disp(2)),N);
    tayPI{j} = tPr;
    PI{j} = CoeffPr;
    Average(j) = mean(PI{j});

    clear Z1 t1
end
%% Построение бифуркационной диаграммы

figure;
subplot(211)
hold on; box on; grid on;
for k = 1:le_vect_ex
    marker = '.k';
    if Average(k)<0
        marker = '.r';
    end 
    plot(vect_ex(k),MaxKSIx{k},marker,'MarkerSize',18);
end
    ff = gca; 
    ff.FontName = 'Times New Roman';
    ff.FontSize = 20;
    xlabel('\epsilon');
    ylabel('EXTR[\xi_{\it x}]')
    ylim([-0.03,0.03]);
    xlim('padded')
subplot(212)
hold on; box on; grid on;
plot(vect_ex,Average,'.-r','MarkerSize',18,'Color','#7E2F8E')
yline(0,'--k','LineWidth',2)
    ff = gca; 
    ff.FontName = 'Times New Roman';
    ff.FontSize = 20;
    xlabel('\epsilon');
    ylabel('Avg(\Lambda)')
    ylim([-1,1.5]);
    xlim('padded')
    
%% Совместные графики
figure; 
hold on; box on; grid on;
plot(vect_ex,Average1,'.-','MarkerSize',22,'LineWidth',2,'Color','#0072BD')
plot(vect_ex,Average2,'.-','MarkerSize',22,'LineWidth',2,'Color','#D95319')
plot(vect_ex,Average3,'.-','MarkerSize',22,'LineWidth',2,'Color','#7E2F8E')
yline(0,'--k','LineWidth',2)
    ff = gca;
    ff.FontName = 'Times New Roman';
    ff.FontSize = 20;
    xlabel('\epsilon');
    ylabel('Avg(\Lambda)')   
    ylim([-1,1.5]);
    xlim('padded')
    legend('Одна опора','Две опоры','Три опоры')
%% Теперь можно смотреть различные реализации 

%EigVal = Calc_EigenValues(A0,A1,A2);
%% Графики перемещений и углов поворота
    tNS = 30; % Количство последних оборотов
    ind = 24;
Moving_the_center(time{ind},ZZ{ind},TN,tNS,disp,angle,N)

%% Траектория центра диска со стробоскопическим отображением
    tNS = 300; % Количство последних оборотов
    ind = 24;
Phase_trajectory_2D(time{ind},ZZ{ind},TN,tNS,disp,N)

%% Квазифазовая траектория центра диска
    tNS = 300; % Количство последних оборотов
    ind = 35;
Quasi_phase_trajectorys(time{ind},ZZ{ind},TN,tNS,disp)

%% Спектральный анализ перемещений
    tNS = 300; % Количество последних оборотов
    ind = 35;
Spectral_analysis(time{ind},ZZ{ind},tNS,TN,disp)

%% Реакции в упругодемпфирующих опорах
    tNS = 20; % Количство последних оборотов
    ind = 24;
Reactions(time{ind},ZZ{ind},TN,tNS,NumbOp,disp,disp_vel,kappa_j,h_j,d_j,f_j,N,eD)
%% Коэффициент прецессии
    tNS = 200; % Количство последних оборотов
    ind = 24;
Precession_coefficient(time{ind},ZZ{ind},TN,tNS,disp,N)

%% Построение положения узлов в 3D
    tNS = 50; % Количство последних оборотов
    ind = 24;
plot_Disp_3D(time{ind},ZZ{ind},TN,tNS,m,index_disp,jC,System_name)



