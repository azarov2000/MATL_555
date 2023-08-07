close all
clc
clear
%% Дано:
N = 14.5;             % Скорость вращения ротора
m = 6;              % Количество участков
NumbOp=60:120:300;  % Расположение опор (degrees)

zDisk=0.7;      % Кооридната расположения диска
System_name=4;  % Тип системы Заделка – заделка(4), консоль (3), Заделка – шарнир (1)
tNS=300;        % Количество оборотов для фазовой трактории
tNS2=300;       % Количество оборотов для стробоскопического отображения
%--------------------------------------------------------------------------

% Размерные параметры
ro = 7800;              % [кг/м^3] - плотность материала
Elastic = 2.1*10^11;    % [Па] - модуль упругости материала 
l = 0.7;                % [м] - длина стержня
d = 20*10^-3;           % [м] - диаметр стержня 
D = 0.2;                % [м] - диаметр диска
a = 10*10^-3;           % [м] - толщина диска
M = ro*(pi*D^2/4)*a;    % [кг] - масса диска
B = M*D^2/16;           % [кг*м^2] - физический момент инерции диска
IR = pi*d^4/64;         % [м^4] - геометрический момент инерции сечения стержня
mR = ro*pi*d^2/4;       % [кг/м] - погонная масса стержня

beta = B/(M*l^2);       % безразмерный коэффициент момента инерции диска
betaR = ro*IR/(M*l);    % безразмерный коэффициент поворотной инерции 
muR = mR*l/M;           % безразмерный коэффициент массы стержня

zi=0.005;               % коэффициент внутреннего линейного демпфирования ротора
ze=0.025;               % коэффициент внешнего линейного демпфирования ротора
Ze=0.0025;              % коэффициент внешнего линейного демфирования диска 
Zte=0.05;               % коэффициент внешнего углового демфирования диска

kappa_s = 850;          % коэффициент жесткости ограничителей
D_s = 0.06;             % коэффициент демпфировния ограничителей   
f_s = 0.1;              % коэффициент сухого трения ограничителей 
eta_s = (10*10^-3)/l;  % коэффициент зазора
nn = 1;                 % cтепень демпфирования
eD = D/l;               % коэффициент диаметра диска
e =(1*10^-3)/l;         % коэффициент эксцентриситета
%gd = 9.81*M*l^2/(Elastic*IR); % коэффициент силы тяжести
gd = 0;
TT = sqrt((M*l^3)/(Elastic*IR));
%%
h=1/m; % Граничные случаи задание постановки диска
s=round(zDisk/h);
if s==0
    jC=1;
elseif s==m
    jC=m;
else
    jC=s;
end
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

aksi=kron(I,E+2*zi*N*S);
dksi=kron(Z,E);
ateta=kron(Z,E);
dteta=kron(I,S-2*zi*N*E);


bksi=kron(2*zi*I+2*ze*h*G00Int+2*Ze*G00*EC,E);
eksi=kron(2*betaR*N*(-h*G01Int+G00*DE)+2*beta*N*G01*EC,E)+kron(2*Zte*G01*EC,S);
bteta=kron(2*ze*h*G10Int+2*Ze*G10*EC,E);
eteta=kron(2*betaR*N*(-h*G11Int+G10*DE)+2*beta*N*G11*EC,E)+kron(2*zi*I+2*Zte*G11*EC,S);


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
x0=0.001*[x0;zeros(2*(mnew+1),1);zeros(4*(mnew+1),1)]; % вектор начальных условий

DD=[kron(G00*EC,E) , zeros(2*(mnew+1),2*(mnew+1))   ;   zeros(2*(mnew+1),2*(mnew+1))   ,   kron(G10*EC,E)]; % Матрица коэффициентов [D] при векторе P тильда (правая часть)

 %% Изменение нумирации
%ii=6;% Для m=6; Zс=0.7; (ii=0), если диск расположен посередине
ii = 2; %(для Zc = 0.7 – заделка заделка)
%ii = 0;
%% Решение системы дифферренциальных уравнений
Ind = 29;
tic;
[t,Z0] = ode23t(@(t,Z0) rGap_Multiple_Supports(t,Z0,DD,A0,A1,A2,N,eD,mnew,eta_s,kappa_s,nn,D_s,e,gd,f_s,NumbOp,ii,System_name),T,x0,opt);
toc;
tN=t/TN; % переход от вектора безразмерного времени к вектору количества оборотов
%% Стробоскопическое отображение
tPoincare=(t(1):TN:t(end))'; % вектор времени, разбитый на равные интервалы, равные времени одного оборота
ZxPoincare=interp1(t,Z0(:,mnew+ii+1),tPoincare,'spline'); % Вектор перемещений по oX через каждый оборот 
ZyPoincare=interp1(t,Z0(:,mnew+ii+2),tPoincare,'spline'); % Вектор перемещений по oY через каждый оборот
tPoincareIndex=tPoincare/TN;                              % Вектор количеств оборотов       
II=find(tPoincareIndex>(tPoincareIndex(end)-tNS2));       % Индексы вектора последних оборотов
%
I=find(tN>(tN(end)-tNS)); % Индексы вектора последних оборотов

%% Коэффициент скорости прецессии
timePI = TN*tN(I);
[tPr,CoeffPr] = PrecessionCoeff(timePI,Z0(I,mnew+ii+1),Z0(I,mnew+ii+2),N);
tPI = tPr;
PI = CoeffPr;
PowerSignal = 0.999;
[fPi,RealizationPi] = spectrum_Fig_Power(tPI,PI,0,1,PowerSignal);
Average = mean(PI);

figure('WindowState','maximized');
     subplot(2,1,1)
     plot(tPI,PI)
     xlabel('\tau','FontName','Times New Roman','FontSize',16)
     ylabel('\Lambda','FontName','Times New Roman','FontSize',16)
     title([num2str(nT-tNS2),' < nT < ',num2str(nT),'; ','N = ', num2str(N)])
     ff = gca; 
     ff.FontName = 'Times New Roman';
     ff.FontSize = 20;
     ylim([0 2]);
     grid on;
     
     subplot(2,1,2)
     h1=stem(fPi,RealizationPi);
     xline(1/TN,'--r',[{'Rotation speed '};'(',num2str(1/TN),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('\itf','FontName','Times New Roman','FontSize',16)
     ylabel('Amplitude','FontName','Times New Roman','FontSize',16)
     title(['Spectrum by \Lambda','; ',num2str(PowerSignal*100),'% signal power'],'FontName','Times New Roman','FontSize',16)
     grid on;
     set(get(h1,'BaseLine'),'LineStyle','-');
     set(h1,'MarkerFaceColor','red');
     ff = gca; 
     ff.FontName = 'Times New Roman';
     ff.FontSize = 20;
     grid on;

%% Графики перемещений и углов поворота
figure('WindowState','maximized');
    subplot(211);hold on;box on
    plot(tN,Z0(:,mnew+ii+1),tN,Z0(:,mnew+ii+2))
    yline(-eta_s,'--k',{'support № 2'},'FontName','Times New Roman','FontSize',14);
    yline(eta_s*0.5,'--k',{'supports № 1,3'},'FontName','Times New Roman','FontSize',14);
    legend('\xi_{ x}','\xi_{ y}','FontName','Times New Roman','FontSize',20)
    
    xlabel(['Displacement, ','{\itt} / (2 \pi / N )'],'FontName','Times New Roman','FontSize',16)
    grid on;
    title(['N = ',num2str(N)]);
    ax2 = gca;
    ax2.FontName = 'Times New Roman';
    ax2.FontSize = 20;
    subplot(212);hold on;box on
    plot(tN,Z0(:,2*(mnew+1)+(mnew)+ii+1),tN,Z0(:,2*(mnew+1)+(mnew)+ii+2))
    legend('\vartheta_{ x}','\vartheta_{ y}','FontName','Times New Roman','FontSize',16)
    xlabel(['Angular movements, ','{\itt} / (2 \pi / N )'],'FontName','Times New Roman','FontSize',16)
    grid on;
    ax2 = gca; 
    ax2.FontName = 'Times New Roman';
    ax2.FontSize = 20;
%% 
%% Фазовая траектрия центра диска
figure3 = figure('WindowState','maximized'); 
    hold on
    plot(-Z0(I,mnew+ii+2),Z0(I,mnew+ii+1),'LineWidth',0.5);
    lenPoin=length(II);
    for j=1:1:lenPoin
        if j<=3
            kColor='o';
            kSize=14;
            p=plot(-ZyPoincare(II(j)),ZxPoincare(II(j)),kColor,'MarkerSize',kSize);
            p.MarkerFaceColor = 'k';
            p.MarkerSize = 8;
        else
            kColor='r.';
            kSize=18;
            plot(-ZyPoincare(II(j)),ZxPoincare(II(j)),kColor,'MarkerSize',kSize);
        end
    end
    axis equal
    xlabel('\xi_{\it y}')
    ylabel('\xi_{\it x}','Rotation',0)
    title([num2str(nT-tNS2),' < nT < ',num2str(nT),'; ','N = ',num2str(N)])
    ax1 = gca;
    ax1.FontName = 'Times New Roman';
    ax1.FontSize = 18;
    grid on; hold on; box on;
    
%% Спектральный анализ
tt=(tN(I)*TN);
Zx=(Z0(I,mnew+ii+1))';
Zy=Z0(I,mnew+ii+2)';
ttInt=linspace(tt(1),tt(end),length(tt));
ZZx=interp1(tt,Zx,ttInt);
ZZy=interp1(tt,Zy,ttInt);
PowerSignal=0.999;
[fx,Xspectr]=spectrum_Fig_Power(ttInt,ZZx,0,1,PowerSignal);
[fy,Yspectr]=spectrum_Fig_Power(ttInt,ZZy,0,1,PowerSignal);
%%
figure4 = figure('WindowState','maximized');
     subplot(2,2,1)
     plot(ttInt,ZZx)
     yline(-eta_s,'--k',{'support № 2'},'LineWidth',3,'FontName','Times New Roman','FontSize',14);
     yline(eta_s*0.5,'--k',{'supports № 1,3'},'LineWidth',3,'FontName','Times New Roman','FontSize',14);
     xlabel('\tau','FontName','Times New Roman','FontSize',16)
     ylabel('\xi_{\itx}','FontName','Times New Roman','FontSize',16)
     %title([num2str(nT-tNS),' < nT < ',num2str(nT),'; ','N = ',num2str(N),'; ','\chi = ', num2str(eta_s)])
     title([num2str(nT-tNS),' < nT < ',num2str(nT),'; ','N = ',num2str(N)])
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;
     ylim([-0.025 0.025]);
     grid on;
     
     subplot(2,2,2)
     h1=stem(fx,Xspectr);
     xline(1/TN,'--r',[{'Rotation speed '};'(',num2str(1/TN),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('\itf','FontName','Times New Roman','FontSize',16)
     ylabel('Amplitude','FontName','Times New Roman','FontSize',16)
     title(['Spectrum by oX','; ',num2str(PowerSignal*100),'% signal power'],'FontName','Times New Roman','FontSize',16)
     grid on;
     set(get(h1,'BaseLine'),'LineStyle','-');
     set(h1,'MarkerFaceColor','red');
     xlim padded
     ylim padded
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;
     grid on;
     
    
     subplot(2,2,3)
     plot(ttInt,ZZy)
     xlabel('\tau','FontName','Times New Roman','FontSize',16)
     ylabel('\xi_{\ity}','FontName','Times New Roman','FontSize',16)
     %title([num2str(nT-tNS),' < nT < ',num2str(nT),'; ','N = ',num2str(N),'; ','\chi = ', num2str(eta_s)])
     title([num2str(nT-tNS),' < nT < ',num2str(nT),'; ','N = ',num2str(N)])
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;
     ylim([-0.025 0.025]);
     grid on;
     grid on;
     
     subplot(2,2,4)
     h2=stem(fy,Yspectr);
     xline(1/TN,'--r',[{'Rotation speed '};'(',num2str(1/TN),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('\itf','FontName','Times New Roman','FontSize',16)
     ylabel('Amplitude','FontName','Times New Roman','FontSize',16)
     title(['Spectrum by oY','; ',num2str(PowerSignal*100),'% signal power'],'FontName','Times New Roman','FontSize',16)
     grid on;
     set(get(h2,'BaseLine'),'LineStyle','-');
     set(h2,'MarkerFaceColor','red');
     xlim padded
     ylim padded
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;
     grid on;
     
     %% Построение квазифазовой траектории
    tEta=t/TN;
    NumerOfRev=find(tEta>(tEta(end)-tNS));
    [textrEtaX,XextrEta]=ext(tEta(NumerOfRev),Z0(NumerOfRev,mnew+ii+1));
    [textrEtaY,YextrEta]=ext(tEta(NumerOfRev),Z0(NumerOfRev,mnew+ii+2));
    
    
    MaxKSIx=XextrEta;
    MaxKSIy=YextrEta;
    le_MAXKSIx = length(MaxKSIx);
    le_MAXKSIy = length(MaxKSIy); 
    vect_x = MaxKSIx(1:(le_MAXKSIx-1));
    vect_y = MaxKSIx(2:le_MAXKSIx);
    
   
    figure;
        grid on; hold on; box on;
        p = plot(vect_x(:),vect_y(:),'.');
        p.MarkerSize = 20;
        p.MarkerFaceColor = 'k';
        x = [-0.025 0.025];
        y = [-0.025,0.025];
        pl = line(x,y);
        pl.Color = 'b';
        pl.LineStyle = '--';
        title(['Квазифазовая траектория за ',num2str(tNS),' последних оборотов'] )
        ff = gca; 
        ff.FontName = 'Times New Roman';
        ff.FontSize = 20;
        xlim([-0.025 0.025]);
        ylim([-0.025 0.025]);
        xticks(-0.025:0.005:0.025);
        yticks(-0.025:0.005:0.025);
        xlabel('\xi_{\it x}^{ \itextr}(N)');
        ylabel('\xi_{\it x}^{ \itextr}(N+1)');
        daspect([1 1 1]) 
   
%% Пояснения
figure;
hold on; box on; grid on; 
    plot(ttInt,ZZx,'LineWidth',2);
    xlabel('\tau')
    ylabel('\xi_{\itx}')
    ff = gca;
    ff.FontName = 'Times New Roman'; 
    ff.FontSize = 20; 
    
   
