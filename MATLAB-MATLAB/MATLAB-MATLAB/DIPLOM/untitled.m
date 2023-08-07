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