close all
clc
clear
%% Дано:
global ze
m=6;                % количество участков
System_name=3;      % случай закрепелния ЗАДЕЛКА – СВОБОДНЫЙ КРАЙ
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
        %xline(0.005,'--black','Внутреннее трение','LineWidth',1.5,'FontName','Times New Roman','FontSize',14);
        plot(zi(:),Ncritical(:),'-r.','MarkerSize',20);
        %plot(zi(:),NcriticalCopy(:),'-b.','MarkerSize',20);
legend('%','Заделка – заделка','Консоль','FontName','Times New Roman','FontSize',20)
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

%% Аппроксимация
x=0; y=0;
for j=1:lenzi
    x(j) = max(IMAG{j});
    y(j) = Ncritical(j);
end
OMEGA=x; NN=y;
X1=[OMEGA(1)/zi(1)    , OMEGA(1);
    OMEGA(end)/zi(end), OMEGA(end)];
X1=[OMEGAMAX(1)/zi(1)    , OMEGAMAX(1);
    OMEGAMAX(end)/zi(end), OMEGAMAX(end)];

NEND=[NN(1);
      NN(end)];
  
a=inv(X1)*NEND;

X2=[OMEGAMIN(1)/zi(1)    , OMEGAMIN(1);
    OMEGAMIN(end)/zi(end), OMEGAMIN(end)];

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

%% Поиск констант a,b и a',b'

% Ncrit=a*ze/zi+b
% omega=a'*ze/zi+b'

% Ncritical - вектор критической скорости 99x1
% zi – вектор внутреннего трения

etai = zi';
Ncr = Ncritical';
omega_plus = OMEGAMAX';
omega_minus = OMEGAMIN';
% Матрица проекта
Nmatrix = [ze./etai ones(size(etai))];
OMEGAmatrix = [ze./etai ones(size(etai))];

% Вычисление коэффициентов
koeffN = Nmatrix\Ncr;
 a = koeffN(1);
 b = koeffN(2);
 
koeffOMEGA_plus = OMEGAmatrix\omega_plus;
 aSh_plus = koeffOMEGA_plus(1);
 bSh_plus = koeffOMEGA_plus(2);
 
koeffOMEGA_minus = OMEGAmatrix\omega_minus;
 aSh_minus = koeffOMEGA_minus(1);
 bSh_minus = koeffOMEGA_minus(2);
 
 figure;
 grid on; hold on; box on;
 plot(etai,a*ze./etai+b)
 plot(etai,Ncr,'r.','MarkerSize',20)
 text(0.05,10,['a = ',num2str(a),'; ','b = ',num2str(b)])
 xlabel('eta_{ i}')
 ylabel('N_{critical} = a  eta_{ e} / eta_{i} + b')
 title('N_{critical} vs \eta_{ i}')
 ff = gca;
 ff.FontName = 'Times New Roman';
 ff.FontSize = 20; 

 figure;
 grid on; hold on; box on;
 plot(etai,aSh_plus*ze./etai+bSh_plus)
 plot(etai,aSh_minus*ze./etai+bSh_minus)
 plot(etai,omega_plus,'r.','MarkerSize',20)
 plot(etai,omega_minus,'r.','MarkerSize',20)
 text(0.05,1,['a = ',num2str(aSh_plus),'; ','b = ',num2str(bSh_plus)])
 text(0.05,-1,['a = ',num2str(aSh_minus),'; ','b = ',num2str(bSh_minus)])
 xlabel('\eta_{ i}')
 ylabel('\omega = a  eta_{ e} / eta_{i} + b')
 title('\omega  vs \eta_{ i}')
 ff = gca;
 ff.FontName = 'Times New Roman';
 ff.FontSize = 20; 

%% Совметные графики показателя
figure; 
    grid on; hold on; box on;

    plot(etai,omega_plus,'-r.','MarkerSize',20)
    plot(etai,omega_minus,'-r.','MarkerSize',20)

    plot(etai,omega_plusCopy,'-b.','MarkerSize',20)
    %plot(etai,omega_minusCopy,'-b.','MarkerSize',20)
    xlabel('\eta_{ i}')
    ylabel('\omega','Rotation',0)
    title('Зависимость характеристического показателя от внутреннего трения')
    legend('Заделка – заделка','','Консоль','FontName','Times New Roman','FontSize',20)
    xlim([0 0.102])
    ff = gca;
    ff.FontName = 'Times New Roman';
    ff.FontSize = 20;



