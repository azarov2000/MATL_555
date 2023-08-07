close all
clc
clear
%% Дано:
global ze Ze Zte beta betaR muR
global L E S I Z EC E1 EL h DE 
global Nmass LenN

m=6;             % количество участков
h=1/m;           % длина участка
zDisk=1;         % координата диска 
System_name=3;   % случай закрепелния ЗАДЕЛКА – СВОБОДНЫЙ КРАЙ

s=round(zDisk/h); % Граничные случаи постановки диска
if s==0
    jC=1;
elseif s==m
    jC=m;
else
    jC=s;
end

zi=linspace(0,0.1,101);        %V Коэффициент внутреннего линейного демпфирования ротора
ze=0.025;             %V Коэффициент внешнего линейного демпфирования ротора
Ze=.0025;             %V Коэффициент внешнего линейного демфирования диска 
Zte=0.05;             %V Коэффициент внутреннего углового демфирования диска
beta=0.0051020;       %V Коэффициент физического момента инерции диска
betaR=0.0000357143;   %V Безразмерный коэффициент поворотной инерции 
muR=0.7;              %V Безразмерный коэффициент массы диска
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
Nmax=100;  % Максимальная скорость
Nkol=1000; % Количество  

Nmass=linspace(Nmin,Nmax,Nkol); % Вектор скоростей
LenN=length(Nmass);               % Длина вектора скоростей
%--------------------------------------------------
%% 
lenzi=length(zi);
index=zeros(1,lenzi);
for i=1:lenzi
[index(i)]=NrkVSzi(System_name,zi(i));
end
lenIndex=length(index);
%% Построение графика зависимости критической скорости от величины зазора

ziInterp=0:0.001:0.1;
NInterp=interp1(zi(2:lenIndex),Nmass(index(2:lenIndex)),ziInterp);
%% Nkrit VS etai

figure1 = figure('WindowState','maximized');
axes1 = axes('Parent',figure1);
hold(axes1,'on');grid on; box on;
    for j=2:lenIndex
        if index(j)==0
            index(j)=1;
        end
        plot(zi(j),Nmass(index(j)),'r.','MarkerSize',20)
    end
    plot(ziInterp,NInterp,'-r');

xlabel('\eta_{ i}','FontName','Times New Roman','FontSize',20)
ylabel('N_{ critical} ','FontName','Times New Roman','FontSize',20)
title('Влияние внутреннего трения на критическую скорость','FontName','Times New Roman','FontSize',20)
set(axes1,'FontName','Times New Roman','FontSize',20);

%% Поиск критических Im и Re
for j=2:lenIndex
    N=Nmass(j);U=[];w=[];
    [w] = MatrixOfGreen_Var_N(System_name,Nmass(index(j)),zi(j));
    J=find(abs(w)~=inf);U=w(J);
    WW{j}=U(:);
    ii=find(real(WW{j})>0);
    REAL{j}=real(WW{j}(ii));
    IMAG{j}=imag(WW{j}(ii));
end


%% Nkrit VS omega
figure1 = figure('WindowState','maximized');
axes1 = axes('Parent',figure1);
hold(axes1,'on');grid on; box on;
    for j=3:lenIndex
        if index(j)==0
           index(j)=1;
        end
        plot(IMAG{j},Nmass(index(j)),'r.','MarkerSize',20)
    end

xlabel('\omega','FontName','Times New Roman','FontSize',20)
ylabel('N_{ critical} ','FontName','Times New Roman','FontSize',20)
title('Зависимость критической скорости от показателя','FontName','Times New Roman','FontSize',20)
set(axes1,'FontName','Times New Roman','FontSize',20);

%%
%% Аппроксимация
%%
x=0; y=0;
for j=2:lenIndex
    x(j-1) = (max(IMAG{j}));
    y(j-1) = (Nmass(index((j))));
    zi(j-1)=zi(j);
end
x=x'; y=y';
%%
X=[ones(size(x))  x./zi  x];

a=X\y;

%%


figure1 = figure('WindowState','maximized');
axes1 = axes('Parent',figure1);
hold(axes1,'on');grid on; box on;
    for j=2:lenIndex
        if index(j)==0
           index(j)=1;
        end
        plot((REAL{j}/zi(j)).*IMAG{j},Nmass(index(j)),'r.','MarkerSize',20)
    end

xlabel('(a/\omega+\eta_{ i})','FontName','Times New Roman','FontSize',20)
ylabel('N_{ critical}','FontName','Times New Roman','FontSize',20)
title('Аппроксимация','FontName','Times New Roman','FontSize',20)
set(axes1,'FontName','Times New Roman','FontSize',20);