function [t,Z0,t_for_spectrum,Coeff_for_spectrum,Average]=FINDSOLVE(N,x0,NUMREV)
%% Дано:
m=6;                % Количество участков
NumbOp=60:120:300;  % Расположение опор (degrees)

zDisk=0.7;        % Кооридната расположения диска
System_name=4;  % Тип системы Заделка – заделка(4), консоль (3), Заделка – шарнир (1)
tNS=300;        % Количество оборотов для фазовой трактории
tNS2=300;       % Количество оборотов для стробоскопического отображения
%--------------------------------------------------------------------------

% Размерные параметры
ro = 7800; % [кг/м^3] - плотность материала;
Elastic = 2.1*10^11;    % [Па] - модуль упругости материала; 
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
% eta_s = (10*10^-3)/l;  % коэффициент зазора
eta_s = 0.001;
nn = 1;                 % cтепень демпфирования 
eD = D/l;               % коэффициент диаметра диска
e =(1*10^-3)/l;          % коэффициент эксцентриситета 
gd = 9.81*M*l^2/(Elastic*IR); % коэффициент силы тяжести
%gd = 0;

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

%--------------------------------------------------------------------------

%% Формирование системы дифференциальных уравнений
opt=odeset('AbsTol',1e-7,'RelTol',1e-7); % настройки решателя
TN=2*pi/N;    % время одного периода
nT=3000;      % число периодов - нужно много периодов для получения установившегося решения
tlim=TN*nT;   % полное время исследования
T=[0,tlim];   % интервал полного исследования
DD=[kron(G00*EC,E) , zeros(2*(mnew+1),2*(mnew+1))   ;   zeros(2*(mnew+1),2*(mnew+1))   ,   kron(G10*EC,E)]; % Матрица коэффициентов [D] при векторе P тильда (правая часть)

 %% Изменение нумирации
ii=2;% Для m=6; Zс=0.7; (ii=0), если диск расположен посередине

%% Решение системы дифферренциальных уравнений

[t,Z0]=ode23t(@(t,Z0) rGap_Multiple_Supports(t,Z0,DD,A0,A1,A2,N,eD,mnew,eta_s,kappa_s,nn,D_s,e,gd,f_s,NumbOp,ii,System_name),T,x0,opt);
    tEta=t/TN;
    NumerOfRev=find(tEta>(tEta(end)-NUMREV));
    [textrEtaX,XextrEta]=ext(tEta(NumerOfRev),Z0(NumerOfRev,mnew+ii+1));
    [textrEtaY,YextrEta]=ext(tEta(NumerOfRev),Z0(NumerOfRev,mnew+ii+2));
    MaxKSIx = XextrEta;
    MaxKSIy = YextrEta;
    
    tN = t/TN;
    I = find(tN>(tN(end)-NUMREV)); % Индексы элементов вектора последних оборотов
    tay = TN*tN(I);
    
    % Коэффициент прецессии
    [tPr,CoeffPr] = PrecessionCoeff(tay,Z0(I,mnew+ii+1),Z0(I,mnew+ii+2),N);
    Average = mean(CoeffPr);
    % Спектральный анализ коэффициента прецессии
    PowerSignal = 1;
    step_discretization = 1/50;
    t_for_spectrum = tPr(1):step_discretization:tPr(end);
    t_for_spectrum(end) = [];
    Coeff_for_spectrum = interp1(tPr,CoeffPr,t_for_spectrum);
%     [fPi,RealizationPi,Trand] = spectrum_Fig_Power(t_for_spectrum,Coeff_for_spectrum,0,1,PowerSignal);

end