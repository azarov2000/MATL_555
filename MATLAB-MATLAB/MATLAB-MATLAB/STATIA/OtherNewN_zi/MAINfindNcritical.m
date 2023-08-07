close all
clc
clear
%% Дано:
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
        plot(IMAG{j},Ncritical(j),'r.','MarkerSize',20)
        OMEGAMAX(j)=max(IMAG{j});
        OMEGAMIN(j)=min(IMAG{j});
    end
    plot(OMEGAMIN,Ncritical,'-r',OMEGAMAX,Ncritical,'-r')
xlabel('\omega','FontName','Times New Roman','FontSize',20)
ylabel('N_{ critical} ','FontName','Times New Roman','FontSize',20)
title('Зависимость критической скорости от показателя','FontName','Times New Roman','FontSize',20)
set(axes1,'FontName','Times New Roman','FontSize',20);

%% Аппроксимация
x=0; y=0;
% for j=1:lenzi
%     x(j) = max(IMAG{j});
%     y(j) = Ncritical(j);
% end
% OMEGA=x; NN=y;
% X1=[OMEGA(1)/zi(1)    , OMEGA(1);
%     OMEGA(end)/zi(end), OMEGA(end)];
X1=[OMEGAMAX(1)/zi(1)    , OMEGAMAX(1);
    OMEGAMAX(end)/zi(end), OMEGAMAX(end)];

NEND=[NN(1);
      NN(end)];
a=inv(X1)*NEND;

X2=[OMEGAMIN(1)/zi(1)    , OMEGAMIN(1);
    OMEGAMIN(end)/zi(end), OMEGAMIN(end)];

NEND=[NN(1);
      NN(end)];
b=inv(X2)*NEND;
%%
IntMAX=[OMEGAMAX(1),(OMEGAMAX(1)/zi(1))*a(1)+a(2)*OMEGAMAX(1); OMEGAMAX(end),(OMEGAMAX(end)/zi(end))*a(1)+a(2)*OMEGAMAX(end)];
IntMIN=[OMEGAMIN(1),(OMEGAMIN(1)/zi(1))*b(1)+b(2)*OMEGAMIN(1); OMEGAMIN(end),(OMEGAMIN(end)/zi(end))*b(1)+b(2)*OMEGAMIN(end)];
%%
figure1 = figure('WindowState','maximized');
axes1 = axes('Parent',figure1);
hold(axes1,'on');grid on; box on;
    plot(OMEGAMAX,NN,'*r',OMEGAMIN,NN,'*r')
    plot(IntMAX(:,1),IntMAX(:,2),'-r',IntMIN(:,1),IntMIN(:,2),'-r')
xlabel('\omega','FontName','Times New Roman','FontSize',20)
ylabel('N_{ critical} = \omega(a/\eta_{ i}+b) ','FontName','Times New Roman','FontSize',20)
title('Аппрокисимация','FontName','Times New Roman','FontSize',20)
set(axes1,'FontName','Times New Roman','FontSize',20);

%%
NNNN=(OMEGAMIN./zi)*b(1)+OMEGAMIN*b(2);
figure;
plot(OMEGAMIN,NNNN)
