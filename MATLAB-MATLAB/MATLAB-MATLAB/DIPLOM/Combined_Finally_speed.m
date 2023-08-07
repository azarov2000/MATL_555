function Combined_Finally_speed
clear all
close all

global beta zeta_e zeta_i Zeta_e Zeta_te z_C N nameS

%---------------------------(Дано)-------------------------------------%
beta=0.01;%025;       % коэфф. момента инерции диска относительно оси OZ
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
         col='r.';ms=12;
     else
         col='k.';ms=8;
     end
     plot3(real(W(:,j)),imag(W(:,j)),EF,col,'MarkerSize',ms)
 end
 
 xlabel('Re(\lambda)','Rotation',0)
 ylabel('Im(\lambda)','Rotation',0)
 zlabel('N','Rotation',0)
 title(titlE)
 grid on
 hold off
 %=================================

  %=================================
   % График собственных значений

f2=figure('Name','Im(lambda) VS Re(lambda)');
hold on;box on; 
for j=1:LenN
    if max(real(W(:,j)))>0    
        col='r.';ms=14;
    else
        col='k.';ms=10; 
    end
%     if j==1
%        plot(real(W(:,j)),imag(W(:,j)),'bo','MarkerSize',4,'LineWidth',2);
%     else
        plot(real(W(:,j)),imag(W(:,j)),col,'MarkerSize',ms);
%     end
end
   
    xlabel('RE(\lambda)')
    ylabel('IM(\lambda)','Rotation',0)
    title(titlE)
    %title(titlE,TITLESUP)
    grid on
    hold off
  %=================================

  %=================================
   % Мнимые значения - угловая скорость 
f3=figure('Name','Im(lambda) VS N');
hold on;box on
for j=1:LenN
    if max(real(W(:,j)))>0
        col='r.';ms=14;
    else
        col='k.';ms=10;
    end
  
   plot(NN(j),imag(W(:,j)),col,'MarkerSize',ms);
end
plot(LimN,[0 0],'-k')
plot([Nmin Nmax],[Nmin  Nmax],'-','LineWidth',1,'Color',[0.07 0.62 1])
plot([Nmin Nmax],[Nmin -Nmax],'-','LineWidth',1,'Color',[0.07 0.62 1])
    xlabel('N')
    ylabel('IM(\lambda)','Rotation',0)
    title(titlE)
    %title(titlE,TITLESUP)
    grid on
    hold off
   %=================================
   % Действительные значения - угловая скорость
f3=figure('Name','Re(lambda) VS N');
hold on;box on
for j=1:LenN
    if max(real(W(:,j)))>0
        col='r.';ms=14;
    else
        col='k.';ms=10;
    end
    plot(NN(j),real(W(:,j)),col,'MarkerSize',ms)
end
    xlabel('N')
    ylabel('RE(\lambda)','Rotation',0)
   title(titlE)
    %title(titlE,TITLESUP)
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
     eD=0.4;        % Коэффициент диаметра диска
     e=0.01;       % Коэффициент эксцентриситета 
    gd=0.001488;     %коэффициет силы тяжести
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
opt=odeset('AbsTol',1e-9,'RelTol',1e-9);
TN=2*pi/N;    % время одного периода
nT=2000;      % число периодов
tlim=TN*nT;   % подное время исследования
T=[0,tlim];
x0=[1;1;0;0;0;0;0;0]*0.01; % начальные условия
[t,Y]=ode45(@(t,Y) rGap(t,Y,e,eta_s),T,x0,opt);
tN=t/TN;
TNT=(0:TN:tlim);
[t1,Y1]=ode45(@(t1,Y1) rGap1(t1,Y1,e,eta_s),TNT,x0,opt);
TNN=t1/TN;

%=========================================

TT=linspace(0,tlim,50000);
tnn=TT/TN;
III=find(tnn>(tnn(end)-20));

% Зависимость числа экстремумов от эксцентриситета
% e1=linspace(0,0.04,50);
% le=length(e1);
% LeTNT=length(III);
% YYY=zeros(LeTNT,le);
% Yex=zeros(100,le);
% for ee=1:le
% [te,Ye]=ode45(@(t3,Y3) rGap(t3,Y3,e1(ee),eta_s),TT,x0,opt);
% YYY(:,ee)=Ye(III,1);
% [textr,Yextr]=ext(te(III),Ye(III,1));
% leY(ee)=length(Yextr);
% Yex(:,ee)=[Yextr; zeros(100-leY(ee),1)];
% end
% 
% % Точки Пуанкаре от изменения эксцентриситета
%  f10=figure;        % xi_y VS xu_x       
%     hold on;box on
%     for j=1:le
%         col1='r.';ms1=8;
%      plot(e1(j),Yex(1:leY(j),j),col1,'MarkerSize',ms1)
%     end
%     xlabel('\epsilon')
%     ylabel('\xi_{ x}')
%     axis equal
%     grid on;
% 
% % Зависимость числа экстремумов от зазора
% eta1=linspace(0,0.04,50);
% leta=length(eta1);
% LeTNT=length(III);
% YYYeta=zeros(LeTNT,leta);
% YexEta=zeros(100,leta);
% for ee=1:leta
% [teta,Yeta]=ode45(@(t3,Y3) rGap(t3,Y3,e,eta1(ee)),TT,x0,opt);
% YYYeta(:,ee)=Yeta(III,1);
% [textrEta,YextrEta]=ext(teta(III),Yeta(III,1));
% leYeta(ee)=length(YextrEta);
% YexEta(:,ee)=[YextrEta; zeros(100-leYeta(ee),1)];
% end
% f10=figure;        % xi_y VS xu_x       
%     hold on;box on
%     for j=1:leta
%         col1='r.';ms1=8;
%      plot(eta1(j),YexEta(1:leYeta(j),j),col1,'MarkerSize',ms1)
%     end
%     xlabel('\eta_{ s}')
%     ylabel('\xi_{ x}')
%     axis equal
%     grid on;



%============================================
% Графики перемещиний и глов поворота
figure;         % Y VS t/TN ...
    subplot(211);hold on;box on
    plot(tN,Y(:,1),tN,Y(:,2))
    legend('\xi_{ x}','\xi_{ y}')
    title(title2)
    grid on;
    
    subplot(212);hold on;box on
    plot(tN,Y(:,3),tN,Y(:,4))
    legend('\theta_{ x}','\theta_{ y}')
    xlabel('{\itt} / (2 \pi / N )')
    grid on;


%============================================
%Цикл для записи сил в упругих опорах
ltt=length(t);  C11=zeros(2,ltt);  F11=zeros(2,ltt);
C22=zeros(2,ltt);  F22=zeros(2,ltt);
for k=1:ltt 
      Vrel1=(N*eD/2+Y(k,6))/abs(N*eD/2+Y(k,6));
      Vrel2=(N*eD/2-Y(k,6))/abs(N*eD/2-Y(k,6));
    C11(:,k)=[1;f_s*Vrel1];
    C22(:,k)=[1;f_s*Vrel2];
     fh1(k)=-((Y(k,1)-eta_s)>=0)*(kappa_s*(Y(k,1)-eta_s)^nn+D_s*Y(k,5));
     fh2(k)=((-Y(k,1)-eta_s)>=0)*(kappa_s*(-Y(k,1)-eta_s)^nn+D_s*Y(k,5));
    F11(:,k)=fh1(k)*C11(:,k);
    F22(:,k)=fh2(k)*C22(:,k);
end
tNS=100; % Число рассматриваемых периодов
I=find(tN>(tN(end)-tNS));% tN>1000-50, т.е. берутся элементы 950-1000
II=find(TNN>(TNN(end)-tNS));
%============================================
% График траектории ценра диска
 figure;        % xi_y VS xu_x       
    hold on;box on
    plot(Y(I,2),Y(I,1),Y1(II,2),Y1(II,1),'r.','MarkerSize',12)
    axis equal
    xlabel('\xi_{ y}')
    ylabel('\xi_{ x}','Rotation',0)
   % title('\xi_{ y} VS \xi_{ x}',TITLESUP)
    grid on;
figure;
comet(-Y(I,2),Y(I,1))
%============================================
% График сил от действия упругой опоры 1 на полном промежутке 
figure;
    subplot(211);hold on;box on;grid on
    plot(tN,F11(1,:))
    ylabel('F1_{ s x}')
    title(title2,TITLESUP)
    
    subplot(212);hold on;box on;grid on
    plot(tN,F11(2,:))
    xlabel('{\itt} / (2 \pi / N )')
    ylabel('F1_{ s y}')

% График сил от действия упругой опоры 2 на полном промежутке 
figure;
    subplot(211);hold on;box on;grid on
    plot(tN,F22(1,:))
    ylabel('F2_{ s x}')
    title(title2)
    
    subplot(212);hold on;box on;grid on
    plot(tN,F22(2,:))
    xlabel('{\itt} / (2 \pi / N )')
    ylabel('F2_{ s y}')
%============================================
% График сил от действия упругой опоры 1 на конечном промежутке
 figure;         % F_x VS t/TN
    subplot(211);hold on;box on;grid on
    plot(tN(I),F11(1,I))
    ylabel('F1_{ s x}')
    title(title2)
    
    subplot(212);hold on;box on;grid on
    plot(tN(I),F11(2,I))
    xlabel('{\itt} / (2 \pi / N )')
    ylabel('F1_{ s y}')

% График сил от действия упругой опоры 2 на конечном промежутке
 figure;         % F_x VS t/TN
    subplot(211);hold on;box on;grid on
    plot(tN(I),F22(1,I))
    ylabel('F2_{ s x}')
    title(title2)
    
    subplot(212);hold on;box on;grid on
    plot(tN(I),F22(2,I))
    xlabel('{\itt} / (2 \pi / N )')
    ylabel('F2_{ s y}')

%============================================
% Функция формирования дифферециальных уравнений
function S=rGap(t,Y,e, eta_s)
    
        c=[cos(N*t);sin(N*t)]; % учёт эксцентриситета
        P=-e*N*(N)*c;

        Vrel1=(N*eD/2+Y(6))/abs(N*eD/2+Y(6)); % учет трения об опору
        C1=[1;f_s*Vrel1];
        Vrel2=(N*eD/2-Y(6))/abs(N*eD/2-Y(6));
        C2=[1;f_s*Vrel2];
        fh1=-((Y(1)-eta_s)>=0)*(kappa_s*(Y(1)-eta_s)^nn+D_s*Y(5));
        fh2=((-Y(1)-eta_s)>=0)*(kappa_s*(-(Y(1))-eta_s)^nn+D_s*Y(5));
        F=fh1*C1+fh2*C2;
        GGG=[-gd;0]; % учет силы тяжести
    ff=A2\[G00*(P+F+GGG);-Gr0*SS*(P+F+GGG)];  
    a=Y(1:4);b=Y(5:8); % формирование правых частей
    S=[b;-B1*b-B0*a+ff];
end

function S1=rGap1(t1,Y1,e, eta_s)
    
        c=[cos(N*t1);sin(N*t1)]; % учёт эксцентриситета
        P=-e*N*(N)*c;

        Vrel1=(N*eD/2+Y1(6))/abs(N*eD/2+Y1(6)); % учет трения об опору
        C1=[1;f_s*Vrel1];
        Vrel2=(N*eD/2-Y1(6))/abs(N*eD/2-Y1(6));
        C2=[1;f_s*Vrel2];
        fh1=-((Y1(1)-eta_s)>=0)*(kappa_s*(Y1(1)-eta_s)^nn+D_s*Y1(5));
        fh2=((-Y1(1)-eta_s)>=0)*(kappa_s*(-(Y1(1))-eta_s)^nn+D_s*Y1(5));
        F=fh1*C1+fh2*C2;
        GGG=[-gd;0]; % учет силы тяжести
    ff=A2\[G00*(P+F+GGG);-Gr0*SS*(P+F+GGG)];  
    a=Y1(1:4);b=Y1(5:8); % формирование правых частей
    S1=[b;-B1*b-B0*a+ff];
end
end


