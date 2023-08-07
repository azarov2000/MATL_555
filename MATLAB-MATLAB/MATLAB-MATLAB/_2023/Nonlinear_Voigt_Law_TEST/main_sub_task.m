close all
clc
clear
% Кубический закон Фойхта
% Схема: заделка-заделка

    global I E R zeta_V Mom  
    global zeta_e G00Int beta_R G02Int G03Int
    global mu_R


warning('on')
m = 6;
h = 1/m;
N = 2;
Mom = 5;

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
beta_R = (eps_d^2)/16;  %[-]
mu_R = 1;

%% Получение матриц
[G00Int,G01Int,G02Int,G03Int,G20Int,G21Int,...
          G00,G01,G02,G03,G20,G21,J]=get_matrix(m);

%%
I=eye(m+1);         % единичная матрица
E = eye(2);         % единичная матрица
EC = zeros(m+1); EC(end,end) = 1;
Z=zeros(m+1);       % нулевая матрица
R = [0 -1; 1 0];    

%% Исключаем первый и последний узлы
% [G00Int,G01Int,G02Int,G03Int,G20Int,G21Int,...
%           G00,G01,G02,G03,G20,G21,J,I,EC,Z]=matrix_edging(G00Int,G01Int,G02Int,G03Int,G20Int,G21Int,...
%           G00,G01,G02,G03,G20,G21,J,I,EC,Z,m);
%!!!!!!
data = update_data_structure(I, E, R, zeta_V, Mom, zeta_e, G00Int, beta_R, ...
    G02Int, G03Int, mu_R);
%% Заполнение матриц коэффициентов

A0 = kron(I,(E-2*zeta_V*N*R))+Mom*kron(G03Int,R);

A1 = 2*zeta_V*kron(I,E)+2*beta_R*N*kron(G02Int,R)+2*zeta_e*kron(G00Int,E);

A2 = mu_R*kron(G00Int,E)-beta_R*kron(G02Int,E);


% g = [kron(G00Int,E), zeros(length(kron(G00Int,E)));zeros(2,length(kron(G00Int,E))) ,kron(G30Int(end,:),E)];



%% Скорость момент
%  Mom_v = linspace(-10,10,100);
% % N_v = linspace(0,50,100);
% N_v = linspace(-60,60,100);


% figure;
% hold on; box on;
% ff = gca;
% ff.FontSize = 25;
% xlabel('M')
% ylabel('N')
% 
% temp_data = data;
%     for i=1:length(Mom_v)
%     i
%     for j=1:length(N_v)
%         temp_data.Mom = Mom_v(i);
%         temp_data.N = N_v(j);
%         [marker_style,marker_size,res] = get_max_real_lyambda(temp_data);
%         plot(Mom_v(i),N_v(j),marker_style,'MarkerSize',marker_size)
%         temp_data = data;
%         Eig_val_critical{i,j} = res;
%     end
% end

Mom_v = linspace(-5,5,100000);
data.N = N;
temp_data = data;
    for j=1:length(Mom_v)
        temp_data.Mom = Mom_v(j);
        [DET_v(j)] = Func(temp_data);
        temp_data = data;
    end
figure;
grid on;
plot(Mom_v,DET_v)
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
step_N = 0.5;

Argan_diagram (N_start,step_N);


%% Решение ДУ

opt=odeset('AbsTol',1e-7,'RelTol',1e-7); % настройки точности

T=[0,1000];   % интервал полного исследования
x0 = zeros(2*length(A0),1); % Вектор начальных условий
amplitude = 0;
period = 2;
phase = 1;
tic;
[t,Z0] = ode15s(@(t,Z0) solver(t,Z0,A0,A1,A2,g,h,N,kappa_B, eta_B,amplitude,phase,period),T,x0,opt);
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
ksi_x = [zeros(length(ksi_x(:,1)),1),ksi_x,Z0(:,length(A0)-1)];

flag = 1;
for j=2:2:(length(A2)-2)
    for i=1:length(Z0(:,1))
        ksi_y(i,j) = Z0(i,j) + Z0(i,length(A0))*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_y(:,1:2:end) = [];
ksi_y = [zeros(length(ksi_y(:,1)),1),ksi_y,Z0(:,length(A0))];

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
ksi_x_dot = [zeros(length(ksi_x_dot(:,1)),1),ksi_x_dot,Z0(:,2*length(A0)-1)];

flag = 1;
for j=length(A2)+2:2:2*(length(A2))-1
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
box on; grid on; hold on;
plot(t,ksi_x(:,end))


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
pause(1);
%%
[t_m,Mi] = get_moment(t,ksi_x,ksi_y,ksi_x_dot,ksi_y_dot,m,N);
%%
t_m(end) = [];
Mi(end) = [];
%%
figure;
box on; grid on; hold on
plot(t_m,Mi);
xlabel('time')
ylabel('M')
ff = gca;
ff.FontSize = 16;

