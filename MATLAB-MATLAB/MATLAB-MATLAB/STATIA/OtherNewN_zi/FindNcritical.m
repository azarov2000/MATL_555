function [Ncritical]=FindNcritical(System_name,zi)
%% Дано:
global ze Ze Zte beta betaR muR
global L E S I Z EC E1 EL h DE

m=6;              % количество участков
h=1/m;            % длина участка
zDisk=1;          % координата диска 
 

s=round(zDisk/h); % Граничные случаи постановки диска
if s==0
    jC=1;
elseif s==m
    jC=m;
else
    jC=s;
end

ze=0.025;              %V Коэффициент внешнего линейного демпфирования ротора
Ze=.0025;              %V Коэффициент внешнего линейного демфирования диска 
Zte=0.05;              %V Коэффициент внутреннего углового демфирования диска
beta=0.0051020;        %V Коэффициент физического момента инерции диска
betaR=0.0000357143;    %V Безразмерный коэффициент поворотной инерции 
muR=0.7;               %V Безразмерный коэффициент массы диска
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
Nmax=100; % Максимальная скорость
Nkol=100; % Количество узлов 

Nmass=linspace(Nmin,Nmax,Nkol); % Вектор скоростей
LenN=length(Nmass);               % Длина вектора скоростей
%--------------------------------------------------
%% Вычисление корней характеристического уравнения
WW=[];
for j=1:LenN
    N=Nmass(j);U=[];w=[];
     [w] = MatrixOfGreen_Var_N(System_name,N,zi);
    J=find(abs(w)~=inf);U=w(J);
    WW{j}=U(:);
end
%% Обработка
ReMax=0;
for j=1:1:LenN
    ReMax(j)=max(real(WW{j}));
end
%%
Npoints=linspace(0,100,10000);
ZEROFUN=@(Npoints) interp1(Nmass(:),ReMax(:),Npoints);
Ncritical=fzero(ZEROFUN,20);




