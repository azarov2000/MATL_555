close all
clc
clear
%% Дано:
global zi ze Ze Zte beta betaR muR
global L E S I Z EC E1 EL h DE

m = 6;              % количество участков
h = 1/m;            % длина участка
zDisk = 0.7;        % координата диска 
System_name = 4;    % случай закрепелния

s = round(zDisk/h); % Граничные случаи постановки диска
if s == 0
    jC = 1;
elseif s == m
    jC = m;
else
    jC = s;
end

zi=0.005;            %v Коэффициент внутреннего линейного демпфирования ротора
ze=0.025;            %v Коэффициент внешнего линейного демпфирования ротора
Ze=.0025;            %v Коэффициент внешнего линейного демфирования диска 
Zte=0.05;            %v Коэффициент внешнего углового демфирования диска
beta=0.0051020;      %v Коэффициент момента инерции диска
betaR=0.0000357143;  %v Безразмерный коэффициент поворотной инерции 
muR=0.7;             %v Безразмерный коэффициент массы диска


L=m+1;
E=eye(2); S=[0 1;-1 0];
I=eye(L); Z=zeros(L);
EC=zeros(L);EC(jC+1,jC+1)=1;
E1=zeros(L);E1(1,1)=1;
EL=zeros(L);EL(L,L)=1;
DE=EL-E1;
%--------------------------
%% Задание скоростей вращения

Nmin = 0;   % Минимальная скорость
%Nmax = 14.262495;  % Максимальная скорость
Nmax = 20;  % Максимальная скорость
Nkol = 100; % Количество узлов 

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
%% Построение трехмерной диаграммы Аргана
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
% f3=figure;
% subplot(211);hold on;box on;grid on
% subplot(212);hold on;box on;grid on
% for j=1:LenN
%     N=Nmass(j);
%     if max(real(WW{j}))>0
%         cc='r.';ms=16;
%     else
%         cc='k.';ms=8;
%     end
%             subplot(211)
%             plot(N,imag(WW{j}),cc,'MarkerSize',ms);
%             subplot(212)
%             plot(N,real(WW{j}),cc,'MarkerSize',ms);
% end
% subplot(211)
% xlabel('N')
% ylabel('Im(\lambda)')
% title('N(Im(\lambda))')
% subplot(212)
% ylabel('Re(\lambda)')
% xlabel('N')
% title('N(Im(\lambda))')

%% Зависимость мнимой части от скорости вращения
f4=figure;
hold on; box on; grid on
for j=1:LenN
    for i=1:48
    N=Nmass(j);
    if real(WW{j}(i))>=0
        cc='r.';ms=16;
    else
        cc='k.';ms=8;
    end
            plot(N,real(WW{j}(i)),cc,'MarkerSize',ms);
    end
end
xlabel('N')
ylabel('Im(\lambda)')
title('Im(\lambda) vs N')
