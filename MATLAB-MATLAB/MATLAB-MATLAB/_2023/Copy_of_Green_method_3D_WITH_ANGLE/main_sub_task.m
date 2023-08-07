close all
clc
clear

    global I E zeta_V S G02Int Mom Fg0 Z G12Int Fg1
    global Fg1_ticks kappa_B G00Int zeta_e beta F00 F10 F02 F12
    global eta_B F12_ticks G10Int G30Int G11Int G01Int F10_ticks
    global mu_B P T

m = 6;
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
% beta = (eps_d^2)/16;    %[-]
beta = 0;
% Параметры краевого условия
m_B = 1;                %[кг]
Stiff_beam = (3*Elastic*J)/l^3;
k_B = Stiff_beam/10;              %[N/m];
%d_B = ///              %

mu_B = m_B/(ro*(pi*d^2/4)*l*l);
kappa_B = k_B*l^3/C;
%kappa_B = 10000
%eta_B = d_B*l/(2*sqrt(ro*(pi*d^2/4)*l*C));

eta_B = 0.01;          %[-]

%% Получение матриц

[F00,F02,Fg0,F10,F12,Fg1,...
F10_ticks,F12_ticks,Fg1_ticks,...
T,P,...
G00Int,G10Int,G01Int,G11Int,G02Int,G12Int,G30Int]=get_matrix_new(m);
%%
I=eye(m+1);         % единичная матрица
E = eye(2);         % единичная матрица
Z=zeros(m+1);       % нулевая матрица
S = [0 1; -1 0];    

%% Исключаем первый и послений узлы
Lim = m+1;
F00 = F00(2:Lim);
F02 = F02(2:Lim);
Fg0 = Fg0(2:Lim);
F10 = F10(2:Lim);
F12 = F12(2:Lim);
Fg1 = Fg1(2:Lim);
T = T(2:Lim);
P = P(2:Lim);

G00Int = G00Int(2:Lim,2:Lim);
G10Int = G10Int(2:Lim,2:Lim);
G01Int = G01Int(2:Lim,2:Lim); 
G11Int = G11Int(2:Lim,2:Lim);
G02Int = G02Int(2:Lim,2:Lim);
G12Int = G12Int(2:Lim,2:Lim);
G30Int = G30Int(2:Lim,2:Lim);
I = I(2:Lim, 2:Lim);
Z = Z(2:Lim, 2:Lim);





%%
A0 = [kron(I,E+2*zeta_V*N*S), -Mom*kron(G02Int,E), -12*Mom*kron(Fg0,S);
      kron(Z,E), kron(I,(E+2*zeta_V*N*S))+Mom*kron(G12Int,S), -12*Mom*kron(Fg1,E);
      kron(Z(end,:),E), -Mom*kron(T',E), -12*Mom*Fg1_ticks*S+12*(E+2*zeta_V*N*S)+6*Mom*S+kappa_B*E];



A1 = [2*zeta_V*kron(I,E)+2*zeta_e*kron(G00Int,E), -2*beta*N*kron(G01Int,E), 2*zeta_e*kron(F00,E)-2*beta*N*kron(F02,S);
      -2*zeta_e*kron(G10Int,S), 2*zeta_V*kron(I,E)+2*beta*N*kron(G11Int,S), -2*zeta_e*kron(F10,S)-2*beta*N*kron(F12,E);
      2*zeta_e*kron(G30Int(end,:),E), -2*beta*N*kron(P',E), 2*zeta_e*F10_ticks*E-2*beta*N*F12_ticks*S+24*zeta_V*E+2*eta_B*E];



A2 = [kron(G00Int,E), beta*kron(G01Int,S), kron(F00,E)-beta*kron(F02,E);
      -kron(G10Int,S), beta*kron(G11Int,E), -kron(F10,S)+beta*kron(F12,S);
      kron(G30Int(end,:),E), beta*kron(P',S), mu_B*E+F10_ticks*E-beta*F12_ticks*E];


g_psi = kron(G00Int,E);
g_theta = -kron(G10Int,S);
g_X = kron(G30Int(end,:),S);

%%
%% Поиск собственных значений
EigVal = polyeig(A0,A1,A2);

figure;
hold on; grid on; box on;
ff = gca;
ff.FontSize = 20;
ff.FontName = 'Times New Roman';

flagNonStab = 0;
for i=1:length(EigVal)
    markerColor = '.k';
    markerSize = 16;
    if real(EigVal(i))>0
        markerColor = '.r';
        markerSize = 25;
        flagNonStab = flagNonStab+1;
    end
    plot(real(EigVal(i)),imag(EigVal(i)),markerColor,'MarkerSize',markerSize)
end
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')

if flagNonStab>0
    title(['Система НЕустойчива ','(max(Re \lambda) = ',num2str(max(real(EigVal))),')'])
else
    title(['Система усточива ','(max(Re \lambda) = ',num2str(max(real(EigVal))),')'])
end


%% Поиск критической скорости
N_start = 0;
step_N = 0.01;

Argan_diagram (N_start,step_N);


%% Решение ДУ

opt=odeset('AbsTol',1e-5,'RelTol',1e-5); % настройки точности

T=[0,1000];   % интервал полного исследования
x0 = zeros(2*length(A0),1); % Вектор начальных условий
amplitude = 0;
period = 2;
phase = 1;
tic;
[t,Z0] = ode23t(@(t,Z0) solver(t,Z0,A0,A1,A2,g_psi,g_theta,g_X,h,N,kappa_B,eta_B,amplitude,phase,period),T,x0,opt);
toc;


%% Формирование ksi_x и ksi_y в каждом узле
flag = 1;
for j=1:2:(length(A2)/2-2)
    for i=1:length(Z0(:,1))
        ksi_x(i,j) = Z0(i,j) + Z0(i,length(A0)-1)*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
%%
ksi_x(:,2:2:end) = [];
ksi_x = [zeros(length(ksi_x(:,1)),1),ksi_x,Z0(:,length(A0)-1)];
%%
flag = 1;
for j=2:2:((length(A2)/2+1)-2)
    for i=1:length(Z0(:,1))
        ksi_y(i,j) = Z0(i,j) + Z0(i,length(A0))*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
%%
ksi_y(:,1:2:end) = [];
ksi_y = [zeros(length(ksi_y(:,1)),1),ksi_y,Z0(:,length(A0))];

z_vect = 0:h:1;

%% Формирование скоростей ksi_x и ksi_y в каждом узле
flag = 1;
for j=length(A2)+1:2:2*(length(A2))-length(A2)/2-2
    for i=1:length(Z0(:,1))
        ksi_x_dot(i,j) = Z0(i,j) + Z0(i,2*length(A0)-1)*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_x_dot(:,1:length(A2))= [];
ksi_x_dot(:,2:2:end) = [];
ksi_x_dot = [zeros(length(ksi_x_dot(:,1)),1),ksi_x_dot,Z0(:,2*length(A0)-1)];

flag = 1;
for j=length(A2)+2:2:2*(length(A2))-length(A2)/2-1
    for i=1:length(Z0(:,1))
        ksi_y_dot(i,j) = Z0(i,j) + Z0(i,2*length(A0))*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_y_dot(:,1:length(A2)+1)= [];
ksi_y_dot(:,2:2:end) = [];
ksi_y_dot = [zeros(length(ksi_y_dot(:,1)),1),ksi_y_dot,Z0(:,2*length(A0))];

z_vect = 0:h:1;

%%
figure;
hold on; box on; grid on;
plot3(z_vect,ksi_x(30,:),ksi_y(30,:))

ff = gca;
ff. FontSize = 16; 
xlabel('z')
ylabel('x')
zlabel('y')

%% 
t_t_vect = linspace(t(end)/1.5,t(end),10000);
ksi_0_x = @(t_t) (amplitude*sin((2*pi/period)*t_t+phase));
ksi_0_y = @(t_t) (amplitude*cos((2*pi/period)*t_t+phase));

for i=1:length(t_t_vect)
    ksi_0_x_(i) = ksi_0_x(t_t_vect(i));
    ksi_0_y_(i) = ksi_0_y(t_t_vect(i));
    ksi_X_1(i) = interp1(t,ksi_x(:,end),t_t_vect(i));
    ksi_Y_1(i) = interp1(t,ksi_y(:,end),t_t_vect(i));
end
%%
figure1 = figure;
axes1 = axes('Parent',figure1);
set(axes1,'DataAspectRatio',[1 1 1]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'on');

plot(ksi_0_x_,ksi_0_y_)
plot(ksi_X_1,ksi_Y_1)

xlabel('\xi_{x}')
ylabel('\xi_{y}','Rotation',0)
legend('Support','Beam')
ff = gca;
ff.FontSize = 16;

%% Трехмерный график перемещений
for i=1:length(z_vect)
    Z_3D{i} = kron(ones(length(ksi_x),1),z_vect(i));
end
index_for_time = find(t<=t_t_vect(end)&t>=t_t_vect(1));
figure;
hold on; box on; grid on;
view([45.9000000862685 13.7999998741406])
for i=1:length(z_vect)
    plot3(Z_3D{i}(index_for_time),ksi_x(index_for_time,i),ksi_y(index_for_time,i),'.-k')
end
ff = gca;
ff.FontSize = 30;

%%
figure;
hold on; box on; grid on
subplot(2, 1, 1)
plot(t,ksi_x(:,end))
plot(t,ksi_y(:,end))
subplot(2, 1, 2)
plot(t,ksi_x_dot(:,end))
plot(t,ksi_y_dot(:,end))

%% Трехмерный график 2
% figure;
% hold on; box on; grid on;
% view([34.8750002358158 46.0463697962445])
% for i=1:length(index_for_time)
%     plot3(z_vect,ksi_x(i,:),ksi_y(i,:),'.-b')
%     drawnow 
%     pause(0.1)
% end
% 
% ff = gca;
% ff. FontSize = 16; 
% xlabel('z')
% ylabel('x')
% zlabel('y')
% pause(1);
%%
[t_m,Mi] = get_moment(t,ksi_x,ksi_y,ksi_x_dot,ksi_y_dot,m,N);
Mi = Mi(1:end-1);
t_m = t_m(1:end-1);
Mi = Mi;
%%
figure;
box on; grid on; hold on
plot(t_m,Mi);
xlabel('time')
ylabel('M')
ff = gca; 
ff.FontSize = 16;


%% Среднее на установившемся
index = find(t_m>=100);

M_TORS_average = mean(Mi(index))


figure; 
plot(t_m(index),Mi(index))



