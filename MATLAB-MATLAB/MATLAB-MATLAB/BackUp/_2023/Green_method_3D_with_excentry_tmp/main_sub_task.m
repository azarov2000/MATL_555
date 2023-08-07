close all
clc
clear
global I E S zeta_V Mom G03Int
global Const_Fg IE Const_Fg_ticks kappa_B 
global zeta_e G00Int beta G02Int Const_F0 Const_F2
global G30E Const_T Const_F0_ticks Const_F2_ticks
global eta_B mu_B

m = 12;
h = 1/m;
N = 2;
Mom = 0;

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
[Const_F0,Const_F2,Const_Fg,Const_Fg30,Const_F0_ticks,Const_F2_ticks,Const_Fg_ticks,Const_T,...
          G00Int,G02Int,G03Int,G30Int]=get_matrix(m);


%%
I=eye(m+1);         % единичная матрица
E = eye(2);         % единичная матрица
Z=zeros(m+1);       % нулевая матрица
S = [0 1; -1 0];    

%% Исключаем первый и послений узлы
Const_F0 = Const_F0(2:end);
Const_F2 = Const_F2(2:end);
Const_Fg = Const_Fg(2:end);
Const_T = Const_T(2:end);
G00Int = G00Int(2:end,2:end);
G02Int = G02Int(2:end,2:end);
G03Int = G03Int(2:end,2:end);
G30Int = G30Int(2:end,2:end);
I = I(2:end,2:end);
Z = Z(2:end);



%% Основные матрицы
IE = kron(I,E);
G30E = kron(G30Int,E);
%%

A0 = [kron(I,E+2*zeta_V*N*S)-Mom*kron(G03Int,S),     -12*Mom*kron(Const_Fg,E);
      12*IE(end-1:end,:),                            -12*Mom*Const_Fg_ticks*S+12*(E+2*zeta_V*N*S)+6*Mom*S+kappa_B*E];

A1 = [2*zeta_V*kron(I,E)+2*zeta_e*kron(G00Int,E)-2*N*beta*kron(G02Int,S),   2*zeta_e*kron(Const_F0,E)-2*beta*N*kron(Const_F2,S);
      2*zeta_e*G30E(end-1:end,:)-2*N*beta*kron(Const_T',S),                 2*zeta_e*Const_F0_ticks*E-2*beta*N*Const_F2_ticks*S+24*zeta_V*E+2*eta_B*E];

A2 = [kron(G00Int,E)-beta*kron(G02Int,E),       kron(Const_F0,E)-beta*kron(Const_F2,E);
      G30E(end-1:end,:)-beta*kron(Const_T',E),  Const_F0_ticks*E-beta*Const_F2_ticks*E+mu_B*E];
 

g = [kron(G00Int,E), zeros(length(kron(G00Int,E)),2);zeros(2,length(kron(G00Int,E))) ,kron(Const_Fg30,E)];
  
%% Поиск собственных значений
% EigVal = polyeig(A0,A1,A2);
% 
% figure;
% hold on; grid on; box on;
% ff = gca;
% ff.FontSize = 20;
% ff.FontName = 'Times New Roman';
% 
% flagNonStab = 0;
% for i=1:length(EigVal)
%     markerColor = '.k';
%     markerSize = 16;
%     if real(EigVal(i))>0
%         markerColor = '.r';
%         markerSize = 25;
%         flagNonStab = flagNonStab+1;
%     end
%     plot(real(EigVal(i)),imag(EigVal(i)),markerColor,'MarkerSize',markerSize)
% end
% xlabel('Re(\lambda)')
% ylabel('Im(\lambda)')
% 
% if flagNonStab>0
%     title(['Система НЕустойчива ','(max(Re \lambda) = ',num2str(max(real(EigVal))),')'])
% else
%     title(['Система усточива ','(max(Re \lambda) = ',num2str(max(real(EigVal))),')'])
% end


%% Поиск критической скорости
N_start = 0;
step_N = 1;

%Argan_diagram (N_start,step_N);


%% Решение ДУ

opt=odeset('AbsTol',1e-5,'RelTol',1e-5,'MassSingular', 'maybe'); % настройки точности

T=[0,1000];   % интервал полного исследования
x0 = zeros(2*length(A0),1); % Вектор начальных условий
amplitude = 0;
period = 2;
phase = 1;
tic;
[t,Z0] = ode23t(@(t,Z0) solver(t,Z0,A0,A1,A2,g,h,N,kappa_B, eta_B,amplitude,phase,period),T,x0,opt);
toc;


%% Формирование ksi_x и ksi_y в каждом узле
flag = 1;
for j=1:2:(length(A2)-2)
    for i=1:length(Z0(:,1))
        ksi_x(i,j) = Z0(i,j) + Z0(i,length(A0)-1)*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_x(:,2:2:end) = [];
ksi_x = [zeros(length(ksi_x(:,1)),1),ksi_x];

flag = 1;
for j=2:2:(length(A2)-2)
    for i=1:length(Z0(:,1))
        ksi_y(i,j) = Z0(i,j) + Z0(i,length(A0))*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_y(:,1:2:end) = [];
ksi_y = [zeros(length(ksi_y(:,1)),1),ksi_y];

z_vect = 0:h:1;

%% Формирование скоростей ksi_x и ksi_y в каждом узле
flag = 1;
for j=length(A2)+1:2:2*(length(A2))-2
    for i=1:length(Z0(:,1))
        ksi_x_dot(i,j) = Z0(i,j) + Z0(i,2*length(A0)-1)*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_x_dot(:,1:length(A2))= [];
ksi_x_dot(:,2:2:end) = [];
ksi_x_dot = [zeros(length(ksi_x_dot(:,1)),1),ksi_x_dot];

flag = 1;
for j=length(A2)+2:2:2*(length(A2))-1
    for i=1:length(Z0(:,1))
        ksi_y_dot(i,j) = Z0(i,j) + Z0(i,2*length(A0))*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_y_dot(:,1:length(A2)+1)= [];
ksi_y_dot(:,2:2:end) = [];
ksi_y_dot = [zeros(length(ksi_y_dot(:,1)),1),ksi_y_dot];

z_vect = 0:h:1;


%% Сравнение PSI
figure;
hold on; box on; grid on;
plot(t,Z0(:,length(A0)-9))
plot(t,Z0(:,length(A0)-7))
plot(t,Z0(:,length(A0)-5))
plot(t,Z0(:,length(A0)-3))
legend('1','2','3','end')


%%
figure;
 hold on; box on; grid on;
 plot3(z_vect,ksi_x(int8(2*length(Z0(:,1))/3),:),ksi_y(int8(2*length(Z0(:,1))/3),:))
 
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
%%
%[t_m,Mi] = get_moment(t,ksi_x,ksi_y,ksi_x_dot,ksi_y_dot,m,N);


