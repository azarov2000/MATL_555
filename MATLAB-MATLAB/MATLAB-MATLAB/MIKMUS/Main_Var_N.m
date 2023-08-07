close all
clc
clear
%% Дано:
global zi ze Ze Zte beta betaR muR
global L E S I Z EC E1 EL h DE

m=6;              % количество участков
h=1/m;            % длина участка
zDisk=0.7;        % координата диска 
System_name=4;    % случай закрепелния
 

s=round(zDisk/h); % Граничные случаи постановки диска
if s==0
    jC=1;
elseif s==m
    jC=m;
else
    jC=s;
end

% Размерные параметры
ro = 7800; % [кг/м^3] - плотность материала;
Elastic = 2.1*10^11; % [Па] - модуль упругости материала; 
l = 0.7;                % [м] - длина стержня
d = 20*10^-3;           % [м] - диаметр стержня 
D = 0.2;                % [м] - диаметр диска
a = 10*10^-3;           % [м] - толщина диска
M = ro*(pi*D^2/4)*a;    % [кг] - масса диска
B = M*D^2/16;           % [кг*м^2] - физический момент инерции диска
IR = pi*d^4/64;         % [м^4] - геометрический момент инерции сечения стержня
mR = ro*pi*d^2/4;       % [кг/м] - погонная масса стержня

beta = B/(M*l^2);       % безразмерный коэффициент момента инерции диска
betaR = ro*IR/(M*l);    % безразмерный коэффициент поворотной инерции 
muR = mR*l/M;           % безразмерный коэффициент массы стержня

zi=0.005;               % коэффициент внутреннего линейного демпфирования ротора
ze=0.025;               % коэффициент внешнего линейного демпфирования ротора
Ze=0.0025;              % коэффициент внешнего линейного демфирования диска 
Zte=0.05;               % коэффициент внешнего углового демфирования диска
%% Задание матриц
L=m+1;
E=eye(2); S=[0 1;-1 0];
I=eye(L); Z=zeros(L);
EC=zeros(L);EC(jC+1,jC+1)=1;
E1=zeros(L);E1(1,1)=1;
EL=zeros(L);EL(L,L)=1;
DE=EL-E1;
%--------------------------
%% Задание скоростей вращения

Nmin=0;   % Минимальная скорость
Nmax=20; % Максимальная скорость
Nkol=602; % Количество узлов 

Nmass=linspace(Nmin,Nmax,Nkol); % Вектор скоростей
LenN=length(Nmass);               % Длина вектора скоростей
%--------------------------------------------------
%% Вычисление корней характеристического уравнения
WW=[];
for j=1:LenN
    N=Nmass(j);U=[];w=[];
     [w] = MatrixOfGreen_Var_N(System_name,N);
    J=find(abs(w)~=inf);U=w(J);
    WW{j}=U(:);
end
%% Построение диаграммы Аргана
f1=figure;
    hold on;box on; grid on
    for j=1:LenN
            if max(real(WW{j}))>=0
                cc='r.';ms=16;
            else 
                cc='k.';ms=8;
            end
            plot(real(WW{j}),imag(WW{j}),cc,'MarkerSize',ms)
    end  
    xlabel('Re( \lambda )')
    ylabel('Im( \lambda )')
    zlabel('N')
%% Построение трехмерной диаграммы Аргна
f2=figure;
    hold on;box on; grid on
    for j=1:LenN
            if max(real(WW{j}))>=0
                cc='r.';ms=16;
            else 
                cc='k.';ms=8;
            end
            plot3(real(WW{j}),imag(WW{j}),Nmass(j)*ones(length(WW{j}),1),cc,'MarkerSize',ms)
    end  
    xlabel('Re( \lambda )')
    ylabel('Im( \lambda )')
    zlabel('N')
%% Построение графиков Re(lambda) vs N и Im(lambda) vs N
f3=figure;
MinRe=-.1;
MaxIm=250;
subplot(211);hold on;box on;grid on
subplot(212);hold on;box on;grid on
for j=1:LenN
    N=Nmass(j);
    if max(real(WW{j}))>0
        cc='r.';ms=16;
    else
        cc='k.';ms=8;
    end
            subplot(211)
            plot(N,imag(WW{j}),cc,'MarkerSize',ms);
            subplot(212)
            plot(N,real(WW{j}),cc,'MarkerSize',ms);
end
subplot(211)
ylabel('Im(\lambda)')
subplot(212)
ylabel('Re(\lambda)')
xlabel('N')

%% Обработка
ReMax=0;
for j=1:1:LenN
    ReMax(j)=max(real(WW{j}));
end

figure;
plot(Nmass(:),ReMax(:))
%%
%FUN = @(Nmass) ReMax;
%%
Npoints=linspace(0,50,0.001);
ZEROFUN=@(Npoints) interp1(Nmass(:),ReMax(:),Npoints);
xx=fzero(ZEROFUN,15)




