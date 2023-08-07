%%
close all
clc
clear
%% Дано:
tic
%--------------------------------------------------------------------------
N=15;               % Скорость вращения ротора
m=6;                % Количество участков
NumbOp = 180:180:180;
%NumbOp = 180:180:180;  % Расположение опор (degrees)
%NumbOp = 36:72:324;

zDisk=1;            % Кооридната расположения диска
System_name=3;      % Тип системы Заделка – заделка(4), консоль (3), Заделка – шарнир (1)
tNS=300;            % Количество оборотов для фазовой трактории
tNS2=300;           % Количество оборотов для стробоскопического отображения

% Размерные параметры
ro = 7800; % [кг/м^3] - плотность материала;
Elastic = 2.1*10^11; % [Па] - модуль упругости материала; 
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
eta_s = (7.5*10^-3)/l;  % коэффициент зазора
nn = 1;                 % cтепень демпфирования 
eD = D/l;               % коэффициент диаметра диска
e =(1*10^-3)/l;               % коэффициент эксцентриситета 
gd = 9.81*M*l^2/(Elastic*IR); % коэффициент силы тяжести
TT = sqrt((M*l^3)/(Elastic*IR));
eta_s = eta_s + (gd/3);
 
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
%% Формирование системы дифференциальных уравнений
opt=odeset('AbsTol',1e-6,'RelTol',1e-6); % настройки решателя
TN=2*pi/N;    % время одного периода
nT=6000;      % число периодов - нужно много периодов для получения установившегося решения
tlim=TN*nT;   % полное время исследования
T=[0,tlim];   % интервал полного исследования
x01=zeros((mnew+1),1); x01(jC,1)=1; x0=kron(x01,ones(2,1));
x0=0.01*[x0;zeros(2*(mnew+1),1);zeros(4*(mnew+1),1)]; % вектор начальных условий
DD=[kron(G00*EC,E) , zeros(2*(mnew+1) , 2*(mnew+1)) ; zeros(2*(mnew+1),2*(mnew+1)) , kron(G10*EC,E)]; % Матрица коэффициентов [D] при векторе P тильда (правая часть)

 %% Изменение нумирации
 ii=6;% Для m=6; Zс=0.7; (ii=0), если диск расположен посередине

%% Вектор по варьируемому параметру
vect_disbalance = linspace(0,e,50);
le_vect_disbalance = length(vect_disbalance);

vect_stiff = 600:50:1200;
le_vect_stiff = length(vect_stiff);

tic
%% 
NUMREV=300;
ind = 29;
for j=1:1:le_vect_stiff
    j
    [t,Z0]=ode23t(@(t,Z0) rGap_Multiple_Supports(t,Z0,DD,A0,A1,A2,N,eD,mnew,eta_s,vect_stiff(j),nn,D_s,e,gd,f_s,NumbOp,ii,System_name),T,x0,opt);
    tEta=t/TN;
    NumerOfRev=find(tEta>(tEta(end)-NUMREV));
    [textrEtaX,XextrEta]=ext(tEta(NumerOfRev),Z0(NumerOfRev,mnew+ii));
    [textrEtaY,YextrEta]=ext(tEta(NumerOfRev),Z0(NumerOfRev,mnew+ii+1));
    MaxKSIx{j}=XextrEta;
    MaxKSIy{j}=YextrEta;
    
    tN{j} = t/TN;
    I = find(tN{j}>(tN{j}(end)-NUMREV)); % Индексы элементов вектора последних оборотов
    tay{j} = TN*tN{j}(I);

    
    x0 = Z0(end,:);
    
    % Коэффициент прецессии
    [tPr,CoeffPr] = PrecessionCoeff(tay{j},Z0(I,mnew+ii),Z0(I,mnew+ii+1),N);
    tayPI{j} = tPr;
    PI{j} = CoeffPr;
    Average(j) = mean(PI{j});
    % Спектральный анализ коэффициента прецессии
    PowerSignal = 0.999;
    [fPi,RealizationPi] = spectrum_Fig_Power(tayPI{j},PI{j},0,1,PowerSignal);
    
    fPI{j} = fPi;
    AmplitudePI{j} = RealizationPi;
    % Сохранение реализаций
    if (j==1)
        ZZ{j} = Z0; time{j} = t;
    end
    clear Z0 t tN
end
toc
%% Построение биффуркационной диаграммы
 figure1 = figure('WindowState','maximized');
 hold on;box on; grid on;
  for j=1:1:le_vect_stiff
      if Average(j)>=0
        plot(vect_stiff(j),(MaxKSIx{j})','k.','MarkerSize',18)
      else 
        plot(vect_stiff(j),(MaxKSIx{j})','r.','MarkerSize',18)
      end
  end
 ff = get(gca);
 %plot(ff1.XLim,[0,0],'k-',ff1.XLim,-[eta_s,eta_s],'k-')
 %plot(ff1.XLim,[0,0],'k-',ff1.XLim,-[eta_s,eta_s],'k-') % одна опора
 plot(ff.XLim,[0,0],'k-',ff.XLim,-[eta_s,eta_s],'k-')               % Одна опора 
 
 xlabel('\kappa','FontName','Times New Roman','FontSize',20)
 ylabel('EXTR [\xi_{\itx}]','FontName','Times New Roman','FontSize',20)   
 title([num2str(nT-NUMREV),' < nT < ',num2str(nT),'; ',...                  % Если "+"
        '\kappa = var(',num2str(vect_stiff(1)),',',num2str(vect_stiff(end)),')'])
 xlim([580 1220]);
 ylim([-0.020 0.020]);
 ff = gca;
 ff.FontName = 'Times New Roman'; 
 ff.FontSize = 20; 
  %% Зависимость Average от значения параметра
    figure5 = figure('WindowState','maximized');
    hold on;box on
    plot(vect_stiff(:),Average(:),'-r')
    plot(vect_stiff(:),Average(:),'r.','MarkerSize',18)
    
    xlabel('\kappa')
    str = '$$ \overline{\Lambda} $$';
    ylabel(str,'interpreter','latex')
    ff = gca; 
    xlim([580 1220]);
    ylim([-1.5 -1])
    ff.FontName = 'Times New Roman';
    ff. FontSize = 20;
    grid on;
    %ylim([-1,1])
%% Просмотр реализаций
 IND = 1;     % Индекс реализации
 tN = time{IND}/TN;
 I = find(tN>(tN(end)-tNS)); % Индексы вектора последних оборотов
 %% Графики перемещений и углов поворота
figure2 = figure('WindowState','maximized');
axes2 = axes('Parent',figure2);
    subplot(211);hold on;box on
    plot(tN,ZZ{IND}(:,mnew+ii),tN,ZZ{IND}(:,mnew+ii+1))
    yline(eta_s-(gd/3),'--k',{'\bfsupport № 1'},'LabelHorizontalAlignment','center','FontName','Times New Roman','FontSize',14);
    yline(-eta_s-(gd/3),'--k',{'\bfsupport №2'},'LabelHorizontalAlignment','center','FontName','Times New Roman','FontSize',14);
    %yline(eta_s*0.5,'--k',{'supports № 1,3'},'FontName','Times New Roman','FontSize',14);
    legend('\xi_{ x}','\xi_{ y}','FontName','Times New Roman','FontSize',20)
    xlabel(['Displacement, ','{\itt} / (2 \pi / N )'],'FontName','Times New Roman','FontSize',16)
    grid on;
    title(['\kappa = ' ,num2str(vect_stiff(IND))])
    ff = gca; 
    ff.FontSize = 20;
    ff.FontName = 'Times New Roman';
    xlim([0 6000])
    ylim([-0.02 0.02])
    subplot(212);hold on;box on
    plot(tN,ZZ{IND}(:,2*(mnew+1)+(mnew)+ii),tN,ZZ{IND}(:,2*(mnew+1)+(mnew+1)+ii))
    
    legend('\vartheta_{ x}','\vartheta_{ y}','FontName','Times New Roman','FontSize',16)
    xlabel(['Angular movements, ','{\itt} / (2 \pi / N )'],'FontName','Times New Roman','FontSize',16)
    grid on;
    ff = gca; 
    ff.FontSize = 20;
    ff.FontName = 'Times New Roman';
    xlim([0 6000])
    ylim([-0.025 0.025])
%% Фазовая траектория
Ratio = [-0.018 0.018]*1;
% Стробоскопическое отображение
tPoincare=(0:TN:tlim)';                                               % вектор времени, разбитый на равные интервалы, равные времени одного оборота
ZxPoincare=interp1(time{IND},ZZ{IND}(:,mnew+ii),tPoincare,'spline');   % Вектор перемещений по oX через каждый оборот 
ZyPoincare=interp1(time{IND},ZZ{IND}(:,mnew+ii+1),tPoincare,'spline'); % Вектор перемещений по oY через каждый оборот
tPoincareIndex=tPoincare/TN;                                          % Вектор количеств оборотов                   
II=find(tPoincareIndex>(tPoincareIndex(end)-tNS2));                   % Индексы вектора последних оборотов


figure3 = figure('WindowState','maximized');
 axes3 = axes('Parent',figure3);
    hold on;box on
    plot(-ZZ{IND}(I,mnew+ii+1),ZZ{IND}(I,mnew+ii),'LineWidth',0.5);
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
    xlabel('\xi_{ y}')
    ylabel('\xi_{ x}','Rotation',0)
     title([num2str(nT-tNS2),' < nT < ',num2str(nT),'; ',...                  % Если "+"
        '\kappa = ' ,num2str(vect_stiff(IND))])
    ff = gca;
    ff.FontName = 'Times New Roman';
    ff.FontSize = 20;
    grid on;    
    xlim(Ratio)
    ylim(Ratio)
    daspect([1 1 1])
%% Фазовая траектория во вращающейся системе
alfa = Average(IND)*N*time{IND}(I);
Ksi_x_rot = ZZ{IND}(I,mnew+ii).*cos(alfa)+ZZ{IND}(I,mnew+ii+1).*sin(alfa);
Ksi_y_rot = -ZZ{IND}(I,mnew+ii).*sin(alfa)+ZZ{IND}(I,mnew+ii+1).*cos(alfa);

alfa_Poin = Average(IND)*N*tPoincare(II);
Poincare_x = ZxPoincare(II).*cos(alfa_Poin)+ZyPoincare(II).*sin(alfa_Poin);
Poincare_y = -ZxPoincare(II).*sin(alfa_Poin)+ZyPoincare(II).*cos(alfa_Poin);
figure;
box on; hold on; grid on;
plot(-Ksi_y_rot,Ksi_x_rot,'LineWidth',0.5)
plot(-Poincare_y,Poincare_x,'.r','MarkerSize',18)
ff = gca;
ff.FontName = ' Times New Roman';
ff.FontSize = 20;
title([num2str(nT-tNS2),' < nT < ',num2str(nT),'; ','\kappa = ',num2str(vect_stiff(IND))])
xlabel('\xi_{\ity}^{ вращ}')
ylabel('\xi_{\itx}^{ вращ}','Rotation',0)
xlim(Ratio)
ylim(Ratio)
daspect([1 1 1])
%% Спектральный анализ для одной из реализаций
tt=(tN(I)*TN);
Zx=(ZZ{IND}(I,mnew+ii))';
Zy=ZZ{IND}(I,mnew+ii+1)';
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
     %yline(-eta_s,'--k',{'support № 2'},'LineWidth',3,'FontName','Times New Roman','FontSize',14);
     %yline(eta_s*0.5,'--k',{'supports № 1,3'},'LineWidth',3,'FontName','Times New Roman','FontSize',14);
     xlabel('\tau','FontName','Times New Roman','FontSize',16)
     ylabel('\xi_{\itx}','FontName','Times New Roman','FontSize',16)
     title(['Moves on oX','; ',num2str(nT-tNS),' < nT < ',num2str(nT),';  ','\kappa = ',num2str(vect_stiff(IND))],'FontName','Times New Roman','FontSize',16)
     ff = gca; 
     ff.FontName = ' Times New Roman';
     ff.FontSize = 20;
     grid on;
     ylim([-0.035 0.035])
     
     subplot(2,2,2)
     h1=stem(fx,Xspectr);
     xline(1/TN,'--r',[{'Rotation speed '};'(',num2str(1/TN),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('Normalized Frequency','FontName','Times New Roman','FontSize',16)
     ylabel('Amplitude','FontName','Times New Roman','FontSize',16)
     title(['Spectrum by oX','; ',num2str(PowerSignal*100),'% signal power'],'FontName','Times New Roman','FontSize',16)
     grid on;
     set(get(h1,'BaseLine'),'LineStyle','-');
     set(h1,'MarkerFaceColor','red');
     xlim padded
     ylim padded
     ff = gca; 
     ff.FontName = ' Times New Roman';
     ff.FontSize = 20;
     
    
     subplot(2,2,3)
     plot(ttInt,ZZy)
     xlabel('\tau','FontName','Times New Roman','FontSize',16)
     ylabel('\xi_{\ity}','FontName','Times New Roman','FontSize',16)
     title(['Moves on oY','; ',num2str(nT-tNS),' < nT < ',num2str(nT),';  ','\kappa = ',num2str(vect_stiff(IND))],'FontName','Times New Roman','FontSize',16)
     grid on;
     ff = gca; 
     ff.FontName = ' Times New Roman';
     ff.FontSize = 20;
     ylim([-0.035 0.035])
     
     subplot(2,2,4)
     h2=stem(fy,Yspectr);
     xline(1/TN,'--r',[{'Rotation speed '};'(',num2str(1/TN),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('Normalized Frequency','FontName','Times New Roman','FontSize',16)
     ylabel('Amplitude','FontName','Times New Roman','FontSize',16)
     title(['Spectrum by oY','; ',num2str(PowerSignal*100),'% signal power'],'FontName','Times New Roman','FontSize',16)
     grid on;
     set(get(h2,'BaseLine'),'LineStyle','-');
     set(h2,'MarkerFaceColor','red');
     xlim padded
     ylim padded
     ff = gca; 
     ff.FontName = ' Times New Roman';
     ff.FontSize = 20;
     %% Спектральный анализ для коэффициента прецессии
indexPressigion = IND;
figure4 = figure('WindowState','maximized');
ff=get(gca);
     subplot(2,1,1)
     plot(tayPI{indexPressigion},PI{indexPressigion})
     xlabel('\tau','FontName','Times New Roman','FontSize',16)
     ylabel('\Lambda','FontName','Times New Roman','FontSize',16)
    title([num2str(nT-NUMREV),' < nT < ',num2str(nT),'; ','\kappa = ',num2str(vect_stiff(indexPressigion))])
     ff = gca; 
     ff.FontName = 'Times New Roman';
     ff.FontSize = 20; 

     grid on;
     
     subplot(2,1,2)
     h1=stem(fPI{indexPressigion},AmplitudePI{indexPressigion});
     xline(1/TN,'--r',[{'Rotation speed '};'(',num2str(1/TN),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('Normalized Frequency','FontName','Times New Roman','FontSize',16)
     ylabel('Amplitude','FontName','Times New Roman','FontSize',16)
     title(['Spectrum by \Lambda','; ',num2str(PowerSignal*100),'% signal power'],'FontName','Times New Roman','FontSize',16)
     grid on;
     set(get(h1,'BaseLine'),'LineStyle','-');
     set(h1,'MarkerFaceColor','red');
     ff = gca; 
     ff.FontName = 'Times New Roman';
     ff.FontSize = 20; 
%% Квазифазовая траектория для прецессии
    [extr_t,extr_lambda] = ext(tayPI{indexPressigion},PI{indexPressigion});
     le_extr_lambda = length(extr_lambda); 
     vect_x = extr_lambda(1:(le_extr_lambda-1));
     vect_y = extr_lambda(2:le_extr_lambda);
    
   
figure;
grid on; hold on; box on;
        p = plot(vect_x(:),vect_y(:),'.');
        p.MarkerSize = 20;
        p.MarkerFaceColor = 'k';
        x = [-30,60];
        y = [-30,60];
        pl = line(x,y);
        pl.Color = 'b';
        pl.LineStyle = '--';
        title(['Квазифазовая траектория за ',num2str(tNS),' последних оборотов'] )
        xlim([-30 60]);
        ylim([-30 60]);
        %xticks(-2:0.5:0);
        %yticks(-2:0.5:0);
        ff = gca; 
        ff.FontName = 'Times New Roman';
        ff.FontSize = 20;
        xlabel('\Lambda_{\itextr}(n)');
        ylabel('\Lambda_{\itextr}(n+1)');
        daspect([1 1 1])
     %% Построение квазифазовой траектории перемещений
     le_MAXKSIx = length(MaxKSIx{IND}); 
     vect_x = MaxKSIx{IND}(1:(le_MAXKSIx-1));
     vect_y = MaxKSIx{IND}(2:le_MAXKSIx);
    
   
figure;
grid on; hold on; box on;
        p = plot(vect_x(:),vect_y(:),'.');
        p.MarkerSize = 20;
        p.MarkerFaceColor = 'k';
        x = [-0.1 0.1];
        y = [-0.1,0.1];
        pl = line(x,y);
        pl.Color = 'b';
        pl.LineStyle = '--';
        title(['Квазифазовая траектория за ',num2str(tNS),' последних оборотов'] )
        xlim([-0.035 0.035]);
        ylim([-0.035 0.035]);
        xticks(-0.040:0.005:0.040);
        yticks(-0.040:0.005:0.040);
        ff = gca; 
        ff.FontName = 'Times New Roman';
        ff.FontSize = 20;
        xlabel('\xi_{\itx extr}(n)');
        ylabel('\xi_{\itx extr}(n+1)');
        daspect([1 1 1])


