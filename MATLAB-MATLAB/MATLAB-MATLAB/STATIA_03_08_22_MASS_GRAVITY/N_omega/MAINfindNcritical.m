close all
clc
clear
%% Дано:
global ze
m=6;             % количество участков
System_name=3;   % случай закрепелния ЗАДЕЛКА – СВОБОДНЫЙ КРАЙ
zi=0.002:0.001:0.1; %V Коэффициент внутреннего линейного демпфирования ротора

%--------------------------------------------------
%% 
lenzi=length(zi);
Ncritical=0;
for i=1:lenzi
[Ncritical1]=FindNcritical(System_name,zi(i));
Ncritical(i)=Ncritical1;
end
%% Построение графика зависимости критической скорости от величины зазора
ziInterp=zi(1):0.001:zi(end);
NInterp=interp1(zi,Ncritical,ziInterp);
%% Nkrit VS etai

figure1 = figure('WindowState','maximized');
axes1 = axes('Parent',figure1);
hold(axes1,'on');grid on; box on;
    for j=1:lenzi
        plot(zi(j),Ncritical(j),'r.','MarkerSize',20)
    end
        plot(ziInterp,NInterp,'-r');

xlabel('\eta_{ i}','FontName','Times New Roman','FontSize',20)
ylabel('N_{ critical} ','FontName','Times New Roman','FontSize',20)
title('Влияние внутреннего трения на критическую скорость','FontName','Times New Roman','FontSize',20)
set(axes1,'FontName','Times New Roman','FontSize',20);

 %% Поиск критических Im и Re
 for j=1:lenzi
     U=[];w=[];
     [w] = MatrixOfGreen_Var_N(System_name,Ncritical(j),zi(j));
     J=find(abs(w)~=inf);U=w(J);
     WW{j}=U(:);
     ii=find(real(WW{j}+0.00001)>=0);
     REAL{j}=real(WW{j}(ii));
     IMAG{j}=imag(WW{j}(ii));
 end
% 
 %% Nkrit VS omega
figure1 = figure('WindowState','maximized');
axes1 = axes('Parent',figure1);
hold(axes1,'on');grid on; box on;
    for j=1:lenzi
        plot(Ncritical(j),IMAG{j},'r.','MarkerSize',20)
        OMEGAMAX(j)=max(IMAG{j});
        OMEGAMIN(j)=min(IMAG{j});
    end
    plot(Ncritical,OMEGAMIN,'-r',Ncritical,OMEGAMAX,'-r')
ylabel('Im(\lambda)','FontName','Times New Roman','FontSize',20)
xlabel('N_{ critical} ','FontName','Times New Roman','FontSize',20)
title('Im(\lambda) vs N_{ critical}, при Re(\lambda) = 0','FontName','Times New Roman','FontSize',20)
set(axes1,'FontName','Times New Roman','FontSize',20);

%% %%% Вычисление коэфф. внутреннего трения для заданных скоростей %%%
N_vector = 11:0.5:20;
vect_zi = interp1(Ncritical(:),zi(:),N_vector);

figure;
grid on; hold on; box on;
for j=1:lenzi
        plot(zi(j),Ncritical(j),'r.','MarkerSize',20)
end
        plot(zi(:),Ncritical(:),'-r')
        plot(vect_zi,N_vector,'b.','MarkerSize',20)
%% Проверка записи
le_vect_zi = length(vect_zi);
for j=1:le_vect_zi
     U = [];w = [];
     [w] = MatrixOfGreen_Var_N(System_name,N_vector(j),0.005);
     J = find(abs(w)~=inf);U = w(J);
     WW{j}=U(:);
     ii=find(real(WW{j})>=0);
     REAL_vector{j}=real(WW{j}(ii));
     IMAG_vector{j}=imag(WW{j}(ii));
end
%% Построение графика
figure;
grid on; hold on; box on;
for j = 1:1:le_vect_zi
        if ~isempty(IMAG_vector{j})
    plot(N_vector(j)/(2*pi),IMAG_vector{j}(1),'r.','MarkerSize',16)
        end
end
    
xlabel('N')
ylabel('Im(\lambda)')



 

