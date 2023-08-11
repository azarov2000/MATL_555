close all
clc
clear
%%
% Cubic Feucht's law
% Scheme: sealing-sealing
warning('on')
m = 6;
h = 1/m;
N = 0;

% Rod Parameters
d = 20*10^-3;           %[m]
l = 0.7;                %[m]
ro = 7800;              %[kg/m^3]
Elastic = 2e11;         %[Pa]
mu_R = 1;

zeta_e = 0.025;         %[-] Не менять!
zeta_V = 0.05;         %[-]  Не менять!
% zeta_VV = 0.005;
% alpha = 0.1;
zeta_VV = 0.05;
% alpha = 0.035;
alpha = 0;

N_z = 0;           % Осевая сила
M_z = 0;           % Крутящий момент

Gamma_0 = 7.9;      % Амплитуда силы

%% Get_matrix
[G00Int,G01Int,G02Int,G00,G01,G02]=get_matrix(m);

%%
I=eye(m+1);        
E = eye(2);         
Z=zeros(m+1);       
R = [0 -1; 1 0];    

%% Exclude the first and last nodes
[G00Int,G01Int,G02Int,G00,G01,G02,I,Z] = matrix_edging(G00Int,G01Int,G02Int,G00,G01,G02,I,Z,m);

      
%% Filling in coefficient matrices
M = get_matr_coeff_M(m);

A0 = kron(I,(E-2*zeta_V*N*R)) - N_z*kron(G02Int,E) + M_z*kron(G01Int,R)*M; % ПРОВЕРКА!!!

A1 = 2*zeta_V*kron(I,E)+2*zeta_e*kron(G00Int,E);

A2 = mu_R*kron(G00Int,E);

A2 = (A2+A2')./2;

%% Writing to the structure
[data] = Filling_structure(I,E,R,zeta_V,zeta_e,G00Int,G00,G01Int,G02Int,M,mu_R,M_z,N_z);
%% Search for eigenvalues
[EigVal] = Calc_EigenValues(A0,A1,A2);

%% Find Natural Frequency
Number_freq = 3;
dat.nat_freq = get_nat_freq(data,Number_freq);

%% Поиск критической скорости

N_start = 0;
step_N = 0.5;

Ncritical = get_N_crit (N_start,step_N,data);


%% Solution of the differential equation
T=[0,50];                   % the interval of the full study
% opt=odeset('AbsTol',1e-6,'RelTol',1e-6);
opt = odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'OutputFcn', @(t, y, flag) odeoutput(t, y, flag, T(1), T(end) - T(1)));


x0 = zeros(2*length(A0),1);  % vector of initial conditions
% Node number from 1 to (m-1)
number_node = round((m-1)/2);

% x0(number_node+2) = 0.00005;
% x0(number_node+3) = 0.00005;


tic;
[t,X] = ode23t(@(t,X) solver(t,X,A0,A1,A2,G02Int,G00,E,zeta_VV,alpha,M,N,R,m,Gamma_0,dat.nat_freq),T,x0,opt);
toc;

%% Formation of ksi_x and ksi_y in each node

for i=1:2:length(A0)
    ksi_X{i} = X(:,i);
end
for i=2:2:length(A0)
    ksi_Y{i} = X(:,i);
end
ksi_X(:,2:2:end) = [];
ksi_Y(:,1:2:end) = [];
ksi_X = [zeros(length(ksi_X{1}),1),ksi_X,zeros(length(ksi_X{1}),1)];
ksi_Y = [zeros(length(ksi_Y{1}),1),ksi_Y,zeros(length(ksi_Y{1}),1)];
z_vect = 0:h:1;

for i=1:length(z_vect)
    Z_coord{i} = ones(length(ksi_X{1}),1) * z_vect(i);
end

%% Spectral analysis
dat.tNS = 20;
dat.number_node = round((m-1)/2);
dat.Power = 1;

dat.ksi_x = ksi_X{dat.number_node+1};
dat.ksi_y = ksi_Y{dat.number_node+1};
dat.t = t;
dat.N = N;

Spectral_analysis(dat)


%% 2D movements
% figure;
% 
% subplot(2,1,1)
% hold on; box on; grid on;
% plot(t,ksi_X{number_node+1})
% ff = gca;
% ff.FontSize = 20;
% xlabel('t')
% ylabel('\xi_{ \it x}','Rotation', 0)
% title(['m = ', num2str(m),...
%         ', N = ', num2str(N),', N_{ crit} = ', num2str(Ncritical),...
%         ', \zeta_{ e} = ', num2str(zeta_e),...
%         ', \zeta_{ V} = ', num2str(zeta_V),...
%         ', \zeta_{ VV} = ', num2str(zeta_VV)])
% 
% subplot(2,1,2)
% hold on; box on; grid on;
% plot(t,ksi_Y{number_node+1})
% ff = gca;
% ff.FontSize = 20;
% xlabel('t')
% ylabel('\xi_{ \it y}','Rotation', 0)
% title(['m = ', num2str(m),...
%         ', N = ', num2str(N),', N_{ crit} = ', num2str(Ncritical),...
%         ', \zeta_{ e} = ', num2str(zeta_e),...
%         ', \zeta_{ V} = ', num2str(zeta_V),...
%         ', \zeta_{ VV} = ', num2str(zeta_VV)])


%% Nodes in 3D
Number_end = 20;
ind = find(t>T(end)-Number_end);
figure;
box on; grid on; hold on;

for i=1:length(Z_coord)
plot3(Z_coord{i}(ind),ksi_X{i}(ind),ksi_Y{i}(ind))
ff = gca;
ff.FontSize = 20;

xlabel('\zeta')
ylabel('\it x')
zlabel('\it y')
end

ylim([min(ksi_X{number_node+1}(ind)) max(ksi_X{number_node+1}(ind))]);
zlim([min(ksi_X{number_node+1}(ind)) max(ksi_X{number_node+1}(ind))]);
disp(['Max_disp = ', num2str(max(ksi_X{number_node+1}(ind)))])



