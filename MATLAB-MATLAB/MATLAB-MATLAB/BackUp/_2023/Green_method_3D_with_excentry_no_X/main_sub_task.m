close all
clc
clear
m = 12;
h = 1/m;
N = 2;
Mom = 2;

% Параметры стержня
d = 20*10^-3;           %[m]
l = 0.7;                %[m]
ro = 7800;              %[kg/m^3]
Elastic = 2e11;         %[Pa]
zeta_e = 0.025;         %[-]
zeta_V = 0.005;         %[-]

J = pi*d^4/64;          %[m^4]
C = Elastic*J;
eps_d = d/l;            %[-]
beta = (eps_d^2)/16;    %[-]
% Параметры краевого условия
m_B = 1;                %[кг]
k_B = 1000;              %[N/m];

mu_B = m_B/(ro*(pi*d^2/4)*l*l);
kappa_B = k_B*l^3/C;
%eta_B = d_B*l/(2*sqrt(ro*(pi*d^2/4)*l*C));
eta_B = 0.01;          %[-]

%% Получение матриц
for j=1:m+1     % проходим по строкам матриц Грина
    for k=1:m+1 % проходим по столбцам матриц Грина
        [G00Int(j,k),G02Int(j,k),G03Int(j,k),G30Int(j,k)]=MatrixOfGreenIntegrate(j,k,h);
    end
end


%%
I=eye(m+1);         % единичная матрица
E = eye(2);         % единичная матрица
Z=zeros(m+1);       % нулевая матрица
S = [0 1; -1 0];    

%% Исключаем первый и послений узлы
G00Int = G00Int(2:m,2:m);
G02Int = G02Int(2:m,2:m);
G03Int = G03Int(2:m,2:m);
I = I(2:m,2:m);
Z = Z(2:m);



%% Основные матрицы
%%

A0 = kron(I,(E+2*zeta_V*N*S))-Mom*kron(G03Int,S);

A1 = 2*zeta_V*kron(I,E)-2*beta*N*kron(G02Int,S)+2*zeta_e*kron(G00Int,E);

A2 = kron(G00Int,E)-beta*kron(G02Int,E);
 
G = kron(G00Int,E);
%% Поиск критической скорости
N_start = 0;
step_N = 1;

%Argan_diagram (N_start,step_N);


%% Решение ДУ

opt=odeset('AbsTol',1e-8,'RelTol',1e-8); % настройки точности

T=[0,1000];   % интервал полного исследования
x0 = zeros(2*length(A0),1); % Вектор начальных условий
amplitude = 0;
period = 2;
phase = 1;
tic;
[t,Z0] = ode23t(@(t,Z0) solver(t,Z0,A0,A1,A2,G,h,N),T,x0,opt);
toc;

%% Визуалиция
ksi_x=[];
ksi_y=[];
for i=1:2:length(A2)
    ksi_x = [ksi_x,Z0(:,i)];
end

for i=2:2:length(A2)
    ksi_y = [ksi_y,Z0(:,i)];
end
ksi_x = [zeros(length(ksi_x(:,1)),1),ksi_x,zeros(length(ksi_x(:,1)),1)];
ksi_y = [zeros(length(ksi_y(:,1)),1),ksi_y,zeros(length(ksi_y(:,1)),1)];

z = 0:h:1;
z = kron(z,ones(length(ksi_y(:,1)),1));

index_t = find(t>T(end)/2 & t<= T(end));

figure;
hold on; box on; grid on; 
for i=1:length(ksi_x(1,:))
    plot3(z(index_t,i),ksi_x(index_t,i),ksi_y(index_t,i),'.r','MarkerSize',5)
end





