%%
close all
clc
clear
%% Дано:
tic
%--------------------------------------------------------------------------
N=11.4;               % Скорость вращения ротора
m=6;                % Количество участков
%NumbOp=0:30:330;    % Расположение опор (degrees)
NumbOp=60:120:300;

zDisk=1;      % Кооридната расположения диска
System_name=3;  % Тип системы Заделка – заделка(4), консоль (3), Заделка – шарнир (1)
tNS=300;        % Количество оборотов для фазовой трактории
tNS2=300;       % Количество оборотов для стробоскопического отображения
%--------------------------------------------------------------------------
h=1/m; % Граничные случаи задание постановки диска
s=round(zDisk/h);
if s==0
    jC=1;
elseif s==m
    jC=m;
else
    jC=s;
end
%--------------------------------------------------------------------------
zi=0.005;            %v Коэффициент внутреннего линейного демпфирования ротора
ze=0.025;            %v Коэффициент внешнего линейного демпфирования ротора
Ze=0.0025;           %v Коэффициент внешнего линейного демфирования диска 
Zte=0.05;            %v Коэффициент внешнего углового демфирования диска
beta=0.0051020;      %v Коэффициент момента инерции диска
betaR=0.0000357143;  %v Безразмерный коэффициент поворотной инерции 
muR=0.7;             %v Безразмерный коэффициент массы диска

%%
L=m+1;
E=eye(2); S=[0 1;-1 0];
I=eye(L); Z=zeros(L);
EC=zeros(L);EC(jC+1,jC+1)=1;
E1=zeros(L);E1(1,1)=1;
EL=zeros(L);EL(L,L)=1;
%--------------------------
G00=zeros(L);G01=zeros(L);
G10=zeros(L);G11=zeros(L);
G00Int=zeros(L);G01Int=zeros(L);
G10Int=zeros(L);G11Int=zeros(L);
for j=1:L
    for k=1:L
     [G00(j,k),G01(j,k),G10(j,k),G11(j,k)]=MatrixOfGreen(j,k,System_name,h);
     [G00Int(j,k),G01Int(j,k),G10Int(j,k),G11Int(j,k)]=MatrixOfGreenInt(j,k,System_name,h);
    end
end
%% Окаямление матриц – убираем заведомо однородные степени свободы
mnew=0;
DE=0;
    switch System_name
    
        case 4 % Заделка – заделка
            G00=G00(2:L-1,2:L-1); G01=G01(2:L-1,2:L-1); G10=G10(2:L-1,2:L-1);
            G11=G11(2:L-1,2:L-1); G00Int=G00Int(2:L-1,2:L-1); G01Int=G01Int(2:L-1,2:L-1);
            G10Int=G10Int(2:L-1,2:L-1); G11Int=G11Int(2:L-1,2:L-1);
            I=I(2:L-1,2:L-1); Z=Z(2:L-1,2:L-1); EC=EC(2:L-1,2:L-1);
            E1=E1(2:L-1,2:L-1); EL=EL(2:L-1,2:L-1);
            mnew=m-2;
            DE=EL-E1;

        case 3 % Заделка – свободный край
            G00=G00(2:L,2:L); G01=G01(2:L,2:L); G10=G10(2:L,2:L);
            G11=G11(2:L,2:L); G00Int=G00Int(2:L,2:L); G01Int=G01Int(2:L,2:L);
            G10Int=G10Int(2:L,2:L); G11Int=G11Int(2:L,2:L);
            I=I(2:L,2:L); Z=Z(2:L,2:L); EC=EC(2:L,2:L);
            E1=E1(2:L,2:L); EL=EL(2:L,2:L);
            mnew=m-1;
            DE=EL-E1;

        case 1 % Заделка – шарнир
            G00=G00(2:L,2:L); G01=G01(2:L,2:L); G10=G10(2:L,2:L);
            G11=G11(2:L,2:L); G00Int=G00Int(2:L,2:L); G01Int=G01Int(2:L,2:L);
            G10Int=G10Int(2:L,2:L); G11Int=G11Int(2:L,2:L);
            I=I(2:L,2:L); Z=Z(2:L,2:L); EC=EC(2:L,2:L);
            E1=E1(2:L,2:L); EL=EL(2:L,2:L);
            mnew=m-1;
            DE=EL-E1;
    end

%% Заполнение матриц коэффициентов
%aksi=kron(I,E)+2*zi*h*N*kron(G00Int,S);
%dksi=kron(Z,E);
%ateta=2*zi*h*N*kron(G10Int,S);
%dteta=kron(I,S);
aksi=kron(I,E+2*zi*N*S);
dksi=kron(Z,E);
ateta=kron(Z,E);
dteta=kron(I,S-2*zi*N*E);


% bksi=2*h*(zi+ze)*kron(G00Int,E)+2*Ze*kron(G00*EC,E);
% eksi=-2*betaR*N*h*kron(G01Int,E)+2*betaR*N*kron(G00*EL,E)-...
%     -2*betaR*N*kron(G00*E1,E)-2*beta*N*kron(G01*EC,E)+...
%     +2*Zte*kron(G01*EC,S);
% bteta=2*h*(zi+ze)*kron(G10Int,E)+2*Ze*kron(G10*EC,E);
% eteta=-2*betaR*N*h*kron(G11Int,E)+2*betaR*N*kron(G10*EL,E)-...
%     -2*betaR*N*kron(G10*E1,E)-2*beta*N*kron(G11*EC,E)+...
%     +2*Zte*kron(G11*EC,S);
bksi=kron(2*zi*I+2*ze*h*G00Int+2*Ze*G00*EC,E);
eksi=kron(2*betaR*N*(-h*G01Int+G00*DE)+2*beta*N*G01*EC,E)+kron(2*Zte*G01*EC,S);
bteta=kron(2*ze*h*G10Int+2*Ze*G10*EC,E);
eteta=kron(2*betaR*N*(-h*G11Int+G10*DE)+2*beta*N*G11*EC,E)+kron(2*zi*I+2*Zte*G11*EC,S);



% cksi=muR*h*kron(G00Int,E)+kron(G00*EC,E);
% fksi=betaR*kron(G00*(E1-EL),S)+betaR*h*kron(G01Int,S)+beta*kron(G01*EC,S);
% cteta=muR*h*kron(G10Int,E)+kron(G10*EC,E);
% fteta=betaR*kron(G10*(E1-EL),S)+betaR*h*kron(G11Int,S)+beta*kron(G11*EC,S);
cksi=kron(muR*h*G00Int+G00*EC,E);
fksi=kron(betaR*(h*G01Int-G00*DE)+beta*G01*EC,S);
cteta=kron(muR*h*G10Int+G10*EC,E);
fteta=kron(betaR*(h*G11Int-G10*DE)+beta*G11*EC,S);

A0=[aksi dksi;ateta dteta];

A1=[bksi eksi;bteta eteta];

A2=[cksi fksi; cteta fteta];

% Поиск корней харатеристического уравнения
[w]=polyeig(A0,A1,A2);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Дано: 
 kappa_s=850;    %v Жесткость ограничителя
 D_s=0.06;     %v Коэффициент демпфировния ограничителя
 f_s=0.1;         %v Коэффициент трения об ограничитель
 %eta_s=(7.5*10^-3)/0.7; %v Начальное расстояние до ограничителей
 eta_s = 0.001;
 nn=1;            %v Степень демпфирования 
 eD=0.286;        %v Коэффициент диаметра диска
 e=0.00143;       %v Коэффициент эксцентриситета 
 gd=0;   %v Коэффициент силы тяжести
 
 
 %% Вектора 
vect_disbalance = linspace(0,e,50);
le_vect_disbalance = length(vect_disbalance);


%% Формирование системы дифференциальных уравнений
opt=odeset('AbsTol',1e-7,'RelTol',1e-7); % настройки решателя
TN=2*pi/N;    % время одного периода
nT=6000;      % число периодов - нужно много периодов для получения установившегося решения
tlim=TN*nT;   % полное время исследования
T=[0,tlim];   % интервал полного исследования
x01=zeros((mnew+1),1); x01(jC,1)=1; x0=kron(x01,ones(2,1));
x0=0.01*[x0;zeros(2*(mnew+1),1);zeros(4*(mnew+1),1)]; % вектор начальных условий
DD=[kron(G00*EC,E) , zeros(2*(mnew+1),2*(mnew+1))   ;   zeros(2*(mnew+1),2*(mnew+1))   ,   kron(G10*EC,E)]; % Матрица коэффициентов [D] при векторе P тильда (правая часть)

 %% Изменение нумирации
 ii=6;% Для m=6; Zс=0.7; (ii=0), если диск расположен посередине
%% Решение системы дифферренциальных уравнений
Ind = 29;
tic;
[t,Z0]=ode23t(@(t,Z0) rGap_Multiple_Supports(t,Z0,DD,A0,A1,A2,N,eD,mnew,eta_s,kappa_s,nn,D_s,e,gd,f_s,NumbOp,ii,System_name),T,x0,opt);
toc;


tN=t/TN; % переход от вектора безразмерного времени к вектору количества оборотов
%% Стробоскопическое отображение
tPoincare=(t(1):TN:t(end))'; % вектор времени, разбитый на равные интервалы, равные времени одного оборота
ZxPoincare=interp1(t,Z0(:,mnew+ii),tPoincare,'spline'); % Вектор перемещений по oX через каждый оборот 
ZyPoincare=interp1(t,Z0(:,mnew+ii+1),tPoincare,'spline'); % Вектор перемещений по oY через каждый оборот
tPoincareIndex=tPoincare/TN;                              % Вектор количеств оборотов                   
II=find(tPoincareIndex>(tPoincareIndex(end)-tNS2));       % Индексы вектора последних оборотов
%
I=find(tN>(tN(end)-tNS)); % Индексы вектора последних оборотов
%% Полярная скорость
ALFAspeed=(Z0(I,4*(mnew+1)+(mnew+1)+ii).*Z0(I,mnew+ii) - Z0(I,4*(mnew+1)+(mnew)+ii).*Z0(I,mnew+ii+1))./((Z0(I,mnew+ii)).^2+(Z0(I,mnew+ii+1)).^2);
PI=ALFAspeed/N;
%% Вычисление среднего значения
lenI = length(I); 
INTERPtN = linspace(tN(I(1)),tN(I(end)),lenI+1);
INTERPpi = interp1(tN(I),PI,INTERPtN);

Average = mean (INTERPpi);
SIGN = sign(Average);
%% Графики перемещений и углов поворота
figure2 = figure('WindowState','maximized');
axes2 = axes('Parent',figure2);
    subplot(211);hold on;box on
    plot(tN,Z0(:,mnew+ii),tN,Z0(:,mnew+ii+1))
    yline(-eta_s,'--k',{'support № 2'},'FontName','Times New Roman','FontSize',14);
    yline(eta_s*0.5,'--k',{'supports № 1,3'},'FontName','Times New Roman','FontSize',14);
    legend('\xi_{ x}','\xi_{ y}','FontName','Times New Roman','FontSize',20)
    
    xlabel(['Displacement, ','{\itt} / (2 \pi / N )'],'FontName','Times New Roman','FontSize',16)
    grid on;    
    title(['\eta_{i} = ',num2str(zi),'; ','\eta_{e} = ',num2str(ze),'; ','\eta_{D \xi} = ',num2str(Ze),'; ','\eta_{D\vartheta} = ',num2str(Zte),'; ','\beta = ',num2str(beta),'; ',...
     '\beta = ',num2str(beta),'; ','\beta_{R} = ',num2str(betaR),'; ','\kappa = ',num2str(kappa_s),'; ','D_{s} = ',num2str(D_s),'; ','f = ',num2str(f_s),'; ',...
        '\chi = ',num2str(eta_s),'; ','\epsilon = ',num2str(e),'; ','\gamma = ',num2str(gd),'; ','Average = ',num2str(Average)])
    subplot(212);hold on;box on
    plot(tN,Z0(:,2*(mnew+1)+(mnew)+ii),tN,Z0(:,2*(mnew+1)+(mnew+1)+ii))
    
    legend('\vartheta_{ x}','\vartheta_{ y}','FontName','Times New Roman','FontSize',16)
    xlabel(['Angular movements, ','{\itt} / (2 \pi / N )'],'FontName','Times New Roman','FontSize',16)
    grid on;
%% 
 %% Фазовая траектрия центра диска
%%
figure3 = figure('WindowState','maximized');
 axes3 = axes('Parent',figure3);
    hold on;box on
    plot(-Z0(I,mnew+ii+1),Z0(I,mnew+ii),'LineWidth',2);
    lenPoin=length(II);
    for j=1:1:lenPoin
        if j<=3
            kColor='o';
            kSize=14;
            p=plot(-ZyPoincare(II(j)),ZxPoincare(II(j)),kColor,'MarkerSize',kSize);
            p.MarkerFaceColor = 'k'
            p.MarkerSize = 8;
        else
            kColor='r.';
            kSize=18;
            plot(-ZyPoincare(II(j)),ZxPoincare(II(j)),kColor,'MarkerSize',kSize);
        end
    end
    axis equal
    xlabel('\xi_{ y}')
    ylabel('\xi_{ x}','Rotation',0)
    title(['\eta_{i} = ',num2str(zi),'; ','\eta_{e} = ',num2str(ze),'; ','\eta_{D \xi} = ',num2str(Ze),'; ','\eta_{D\vartheta} = ',num2str(Zte),'; ','\beta = ',num2str(beta),'; ',...
     '\beta = ',num2str(beta),'; ','\beta_{R} = ',num2str(betaR),'; ','\kappa = ',num2str(kappa_s),'; ','D_{s} = ',num2str(D_s),'; ','f = ',num2str(f_s),'; ',...
        '\chi = ',num2str(eta_s),'; ','\epsilon = ',num2str(e),'; ','\gamma = ',num2str(gd),'; ','Average = ',num2str(Average)])
    grid on;
%%
%% Запись сил в упругой опоре
 ltt=length(tN); PItog=zeros(2,ltt); flag=1;
 
 for alfa1=NumbOp
     for k=1:ltt 
        alfa=(alfa1*pi)/180; % переход к радианам
        % Связь перемещений
        uRxs=Z0(k,mnew+ii)*cos(alfa)+Z0(k,mnew+ii+1)*sin(alfa);
        uRys=-Z0(k,mnew+ii)*sin(alfa)+Z0(k,mnew+ii+1)*cos(alfa);
        % Связь скоростей
        duRxs=Z0(k,4*(mnew+1)+(mnew)+ii)*cos(alfa)+Z0(k,4*(mnew+1)+(mnew+1)+ii)*sin(alfa);
        duRys=-Z0(k,4*(mnew+1)+(mnew)+ii)*sin(alfa)+Z0(k,4*(mnew+1)+(mnew+1)+ii)*cos(alfa);
        % Связь сил
        Sn(flag,k)=(kappa_s*(uRxs-eta_s)+D_s*duRxs)*(uRxs-eta_s>=0)*((kappa_s*(uRxs-eta_s)+D_s*duRxs)>=0);
        St(flag,k)=f_s*((N*eD/2+duRys)/abs(N*eD/2+duRys))*Sn(flag,k);
        MatrRot=[cos(alfa), -sin(alfa); sin(alfa),cos(alfa)];
        PItog(:,k)=MatrRot*[-Sn(flag,k);-St(flag,k)];
     end
     flag=flag+1;
 end
%% Comet
figure;
comet(-Z0(I,mnew+ii+1),Z0(I,mnew+ii))
%% Графики сил в упругих опорах  опоре

lenOp=length(NumbOp); % Автоматическая запись легенд
Legend{1}=zeros(1,lenOp);
for j=1:1:lenOp
    
    Legend{j}=['Опора №', num2str(j)];
end

% Графики сил от действия упругих опор полном промежутке 
figure;
    subplot(211);hold on;box on;grid on
       for j=1:1:lenOp
           plot(tN,Sn(j,:))
       end
       legend(Legend)
       title('Нормальные реакции в опорах')
       ylabel('S_ { n\alpha_{j}}')
    
    subplot(212);hold on;box on;grid on
       for j=1:1:lenOp
           plot(tN,St(j,:))
       end
        legend(Legend)
        title('Тангенциальные реакции в опорах')
        xlabel('{\itt} / (2 \pi / N )')
        ylabel('S_ { t\alpha_{j}}')
        
        
% %% Полярная скорость
% ALFAspeed=(Z0(I,4*(mnew+1)+(mnew+1)+ii).*Z0(I,mnew+ii) - Z0(I,4*(mnew+1)+(mnew)+ii).*Z0(I,mnew+ii+1))./((Z0(I,mnew+ii)).^2+(Z0(I,mnew+ii+1)).^2);
% PI=ALFAspeed/N;
% %% Вычисление среднего значения
% lenI = length(I); 
% INTERPtN = linspace(tN(I(1)),tN(I(end)),lenI+1);
% INTERPpi = interp1(tN(I),PI,INTERPtN);
% 
% Average = mean (INTERPpi);
% SIGN = sign(Average);

%% Спектральный анализ
tt=(tN(I)*TN);
Zx=(Z0(I,mnew+ii))';
Zy=Z0(I,mnew+ii+1)';
ttInt=linspace(tt(1),tt(end),length(tt));
ZZx=interp1(tt,Zx,ttInt);
ZZy=interp1(tt,Zy,ttInt);
PowerSignal=0.999;
[fx,Xspectr]=spectrum_Fig_Power(ttInt,ZZx,0,1,PowerSignal);
[fy,Yspectr]=spectrum_Fig_Power(ttInt,ZZy,0,1,PowerSignal);
%%
figure4 = figure('WindowState','maximized');
ff=get(gca);
     subplot(2,2,1)
     plot(ttInt,ZZx)
     yline(-eta_s,'--k',{'support № 2'},'LineWidth',3,'FontName','Times New Roman','FontSize',14);
     yline(eta_s*0.5,'--k',{'supports № 1,3'},'LineWidth',3,'FontName','Times New Roman','FontSize',14);
     xlabel('\tau','FontName','Times New Roman','FontSize',16)
     ylabel('\xi_{_x}','FontName','Times New Roman','FontSize',16)
     title(['Moves on oX','; ',num2str(nT-tNS),' < nT < ',num2str(nT)],'FontName','Times New Roman','FontSize',16)
     grid on;
     
     subplot(2,2,2)
     h1=stem(fx,Xspectr);
     xline(1/TN,'--r',{'Rotation speed'},'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('f','FontName','Times New Roman','FontSize',16)
     ylabel('Amplitude','FontName','Times New Roman','FontSize',16)
     title(['Spectrum by oX','; ',num2str(PowerSignal*100),'% signal power'],'FontName','Times New Roman','FontSize',16)
     grid on;
     set(get(h1,'BaseLine'),'LineStyle','-');
     set(h1,'MarkerFaceColor','red');
     xlim padded
     ylim padded
     
    
     subplot(2,2,3)
     plot(ttInt,ZZy)
     xlabel('\tau','FontName','Times New Roman','FontSize',16)
     ylabel('\xi_{_y}','FontName','Times New Roman','FontSize',16)
     title(['Moves on oY','; ',num2str(nT-tNS),' < nT < ',num2str(nT),';  '],'FontName','Times New Roman','FontSize',16)
     grid on;
     
     subplot(2,2,4)
     h2=stem(fy,Yspectr);
     xline(1/TN,'--r',{'Rotation speed'},'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('f','FontName','Times New Roman','FontSize',16)
     ylabel('Amplitude','FontName','Times New Roman','FontSize',16)
     title(['Spectrum by oY','; ',num2str(PowerSignal*100),'% signal power'],'FontName','Times New Roman','FontSize',16)
     grid on;
     set(get(h2,'BaseLine'),'LineStyle','-');
     set(h2,'MarkerFaceColor','red');
     xlim padded
     ylim padded
