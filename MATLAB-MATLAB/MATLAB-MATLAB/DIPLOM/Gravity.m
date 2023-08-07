function Gravity
clear all
close all

global beta zeta_e zeta_i Zeta_e Zeta_te z_C N nameS

%---------------------------(Дано)-------------------------------------%
beta=0.005;%025;       % коэфф. момента инерции диска относительно оси OZ
zeta_e=0.025;%25;        % коэфф. внешнего демпфирования вала 
zeta_i=0.05;%05;         % коэфф. внутреннего демпфирования вала
Zeta_e=0.025;%025;       % коэфф. внутреннего линейного демпфирования диска 
Zeta_te=0.05;%025;        % коэфф. внутреннего углового демпфирования диска
 %---------------------------------------------------------------------------%

 z_C=input('Введите координату положения диска (0.3 < z_C < 0.7)=');

 ss='Случай закрепления: Заделка-шарнир-"1" | Шарнир-шарнир-"2" | Заделка-свободный край- "3" | Заделка-заделка- "4" ';

 nameS=input(ss);
 switch nameS
     case 1
         SupportName='Заделка-шарнир';
     case 2 
         SupportName='Шарнир-шарнир';
     case 3
         SupportName='Заделка-свободный край';
     case 4
         SupportName='Заделка-заделка';
 end
 %=================================
titlE=['\beta = ',          num2str(beta),...
          ', \zeta_{ i} = ',   num2str(zeta_i),...
          ', \zeta_{ e} = ',   num2str(zeta_e),...
          ', Z_{ e} = ',       num2str(Zeta_e),...
          ', Z_{ \theta e} = ',num2str(Zeta_te),...
          ', \zeta_{ C} = ',   num2str(z_C)];
TITLESUP=[SupportName];
 %=================================
 Nmin=0;
 Nmax=50;
 nN=100;
 NN=linspace(Nmin,Nmax,nN+1); % Записываем вектор угловых скоростей
 LenN=length(NN);

 %=================================
 % Получение собственных векторов и значений при разных угл. скоростях
V=[];W=[]; S=[];

for j=1:LenN
   N=NN(j);
   [A2,A1,A0]=MATRIX(z_C,N);
   [v,w,s]=polyeig(A0,A1,A2); % v-собс.вектора,w-собс.знач.
   V=[V,v(:)]; W=[W,w(:)]; % формир. матр. собств.век. и с.знач.
end
 %=================================
 
 LimN=[Nmin,Nmax]; 

 %=================================
 % Построение графика собсвенных значений от скорости вращения
 f1=figure('Name','Im(lambda) VS Re(lambda) VS N');
 hold on;box on;

 for j=1:LenN
   EF=NN(j)*ones(8,1);
     if max(real(W(:,j)))>0
         col='r.';ms=8;
     else
         col='k.';ms=4;
     end
     plot3(real(W(:,j)),imag(W(:,j)),EF,col,'MarkerSize',ms)
 end
 
 xlabel('Re(\lambda)','Rotation',0)
 ylabel('Im(\lambda)','Rotation',0)
 zlabel('N','Rotation',0)
 title(titlE,TITLESUP)
 grid on
 hold off
 %=================================

  %=================================
   % График собственных значений

f2=figure('Name','Im(lambda) VS Re(lambda)');
hold on;box on; 
for j=1:LenN
    if max(real(W(:,j)))>0    
        col='r.';ms=8;
    else
        col='k.';ms=4; 
    end
    if j==1
       plot(real(W(:,j)),imag(W(:,j)),'bo','MarkerSize',4,'LineWidth',2);
    else
        plot(real(W(:,j)),imag(W(:,j)),col,'MarkerSize',ms);
    end
end
   
    xlabel('RE(\lambda)')
    ylabel('IM(\lambda)','Rotation',0)
    title(titlE,TITLESUP)
    grid on
    hold off
  %=================================

  %=================================
   % Мнимые значения - угловая скорость 
f3=figure('Name','Im(lambda) VS N');
hold on;box on
for j=1:LenN
    if max(real(W(:,j)))>0
        col='r.';ms=8;
    else
        col='k.';ms=4;
    end
  
   plot(NN(j),imag(W(:,j)),col,'MarkerSize',ms);
end
    xlabel('N')
    ylabel('IM(\lambda)','Rotation',0)
    title(titlE,TITLESUP)
    grid on
    hold off
   %=================================
   % Действительные значения - угловая скорость
f3=figure('Name','Re(lambda) VS N');
hold on;box on
for j=1:LenN
    if max(real(W(:,j)))>0
        col='r.';ms=8;
    else
        col='k.';ms=4;
    end
    plot(NN(j),real(W(:,j)),col,'MarkerSize',ms)
end
    xlabel('N')
    ylabel('RE(\lambda)','Rotation',0)
    title(titlE,TITLESUP)
    grid on
    hold off

%***************************************************************
pause(5)


 N=input('N>N_biff=');

 kappa_s=200;        % Жесткость ограничителя
    D_s=0.05;       % Коэффициент демпфировния ограничителя
    f_s=0.3;        % Коэффициент трения об ограничитель
  eta_s=0.03;       % Начальное расстояние до ограничителя
     nn=1;          % Exponent in the law of nonlinearity
     eD=0.2;        % Коэффициент диаметра диска
      e=0.02;       % Коэффициент эксцентриситета 
      gd=0.001;     %коэффициет силы тяжести
     SS=[0 1;-1 0];  ee=eye(2);  oo=zeros(2,2);


 title2=[
        'N = ',              num2str(N),...
        ', \zeta_{ C} = ',   num2str(z_C),...
        ', \epsilon_{ D} = ',num2str(eD),...
        ', \eta_{ s} = ',    num2str(eta_s),...
        ', \kappa_{ s} = ',  num2str(kappa_s),...
        ', D_{ s} = ',       num2str(D_s),...
        ', f_{ s} = ',       num2str(f_s),...
        ', \epsilon = ',     num2str(e),...
        ];

[A2,A1,A0]=MATRIX(z_C,N);
B1=A2\A1;B0=A2\A0;
[G00,G0s,Gr0,Grs]=GRIN(z_C,z_C,nameS);

%============================================
% Настройка решателя
opt=odeset('AbsTol',1e-10,'RelTol',1e-8);
TN=2*pi/N;    % время одного периода
nT=1000;      % число периодов
tlim=TN*nT;   % подное время исследования
T=[0,tlim];
x0=[1;1;0;0;0;0;0;0]*0.02; % начальные условия
[t,Y]=ode45(@(t,Y) rGap(t,Y),T,x0,opt);
tN=t/TN;


%============================================
% Графики перемещиний и глов поворота
figure;         % Y VS t/TN ...
    subplot(211);hold on;box on
    plot(tN,Y(:,1),tN,Y(:,2))
    legend('\xi_{ x}','\xi_{ y}')
    title(title2,TITLESUP)
    
    subplot(212);hold on;box on
    plot(tN,Y(:,3),tN,Y(:,4))
    legend('\theta_{ x}','\theta_{ y}')
    xlabel('{\itt} / (2 \pi / N )')

%============================================
% Цикл для записи сил в упругой опоре
ltt=length(t);  C=zeros(2,ltt);  F=zeros(2,ltt);
for k=1:ltt 
      Vrel=(N*eD/2-Y(k,6))/abs(N*eD/2-Y(k,6));
    C(:,k)=[1;f_s*Vrel];
     fh(k)=-((Y(k,1)-eta_s)>=0)*(kappa_s*(Y(k,1)-eta_s)^nn+D_s*Y(k,5));
    F(:,k)=fh(k)*C(:,k);
end
tNS=50;
 I=find(tN>(tN(end)-tNS));% tN>1000-50, т.е. берутся элементы 950-1000

%============================================
% График траектории ценра диска
 figure;        % xi_y VS xu_x       
    hold on;box on
    plot(Y(I,2),Y(I,1))
    axis equal
   % at=[num2str(tN(I(end))-tNS),' < {\itt} \rm / (2 \pi / N )< ',num2str(tN(I(end)))];
    xlabel('\xi_{ y}','Rotation',0)
    ylabel('\xi_{ x}','Rotation',0)
    title('\xi_{ y} VS \xi_{ x}',TITLESUP)
    grid on
%
%============================================
% График сил от действия упругой опоры на полном промежутке
figure;
    subplot(211);hold on;box on;grid on
    plot(tN,F(1,:))
    ylabel('F_{ s x}')
    title(title2,TITLESUP)
    
    subplot(212);hold on;box on;grid on
    plot(tN,F(2,:))
    xlabel('{\itt} / (2 \pi / N )')
    ylabel('F_{ s y}')

%============================================
% График сил от действия упругой опоры на конечном промежутке
figure;         % F_x VS t/TN
    subplot(211);hold on;box on;grid on
    plot(tN(I),F(1,I))
    ylabel('F_{ s x}')
    title(title2,TITLESUP)
    
    subplot(212);hold on;box on;grid on
    plot(tN(I),F(2,I))
    xlabel('{\itt} / (2 \pi / N )')
    ylabel('F_{ s y}')

%============================================
% Функция формирования дифферециальных уравнений
function S=rGap(t,Y)
    
        c=[cos(N*t);sin(N*t)]; % учёт эксцентриситета
        P=-e*N*(N*ee-2*Zeta_e*SS)*c;

        Vrel=(N*eD/2-Y(6))/abs(N*eD/2-Y(6)); % учет трения об опору
        C=[1;f_s*Vrel];
        fh=-((Y(1)-eta_s)>=0)*(kappa_s*(Y(1)-eta_s)^nn+D_s*Y(5));
        F=fh*C;

        GGG=[-gd;0]; % учет силы тяжести

    ff=A2\[G00*(P+F+GGG);-Gr0*SS*(P+F+GGG)]; % вектор G 
    a=Y(1:4);b=Y(5:8);
    S=[b;-B1*b-B0*a+ff];
end
end


