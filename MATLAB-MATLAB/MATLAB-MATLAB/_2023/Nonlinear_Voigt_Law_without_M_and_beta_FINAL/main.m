close all
clc
clear
% Cubic Feucht's law
% Scheme: sealing-sealing
warning('on')
m = 6;
h = 1/m;
N = 30;

% Rod Parameters
d = 20*10^-3;           %[m]
l = 0.7;                %[m]
ro = 7800;              %[kg/m^3]
Elastic = 2e11;         %[Pa]

zeta_e = 0.025;         %[-]
zeta_V = 0.005;         %[-]
zeta_VV = 0.005;
mu_R = 1;

%% Get_matrix
[G00Int,G02Int,G00,G02]=get_matrix(m);

%%
I=eye(m+1);        
E = eye(2);         
Z=zeros(m+1);       
R = [0 -1; 1 0];    

%% Exclude the first and last nodes
[G00Int,G02Int,G00,G02,I,Z] = matrix_edging(G00Int,G02Int,G00,G02,I,Z,m);

%% Writing to the structure
[data] = Filling_structure(I,E,R,zeta_V,zeta_e,G00Int,G02Int,mu_R);
      
%% Filling in coefficient matrices

A0 = kron(I,(E-2*zeta_V*N*R));

A1 = 2*zeta_V*kron(I,E)+2*zeta_e*kron(G00Int,E);

A2 = mu_R*kron(G00Int,E);


M = get_matr_coeff_M(m);

%% Search for eigenvalues

[EigVal] = Calc_EigenValues(A0,A1,A2);
%% Поиск критической скорости

N_start = 0;
step_N = 0.5;

Ncritical = get_N_crit (N_start,step_N,data);


%% Solution of the differential equation

opt=odeset('AbsTol',1e-6,'RelTol',1e-6);

T=[0,30];                   % the interval of the full study
x0 = zeros(2*length(A0),1); % vector of initial conditions
% Node number from 1 to (m-1)
number_node = 3;

x0(number_node+2) = 0.00005;
x0(number_node+3) = 0.00005;


tic;
[t,X] = ode23t(@(t,X) solver(t,X,A0,A1,A2,G02Int,E,zeta_VV,M,N,R,m),T,x0,opt);
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

%% 2D movements
figure;

subplot(2,1,1)
hold on; box on; grid on;
plot(t,ksi_X{number_node+1})
ff = gca;
ff.FontSize = 20;
xlabel('t')
ylabel('\xi_{ \it x}','Rotation', 0)
title(['m = ', num2str(m),...
        ', N = ', num2str(N),', N_{ crit} = ', num2str(Ncritical),...
        ', \zeta_{ e} = ', num2str(zeta_e),...
        ', \zeta_{ V} = ', num2str(zeta_V),...
        ', \zeta_{ VV} = ', num2str(zeta_VV)])

subplot(2,1,2)
hold on; box on; grid on;
plot(t,ksi_Y{number_node+1})
ff = gca;
ff.FontSize = 20;
xlabel('t')
ylabel('\xi_{ \it y}','Rotation', 0)
title(['m = ', num2str(m),...
        ', N = ', num2str(N),', N_{ crit} = ', num2str(Ncritical),...
        ', \zeta_{ e} = ', num2str(zeta_e),...
        ', \zeta_{ V} = ', num2str(zeta_V),...
        ', \zeta_{ VV} = ', num2str(zeta_VV)])


%% Nodes in 3D
% figure;
% box on; grid on; hold on;
% 
% for i=1:length(Z_coord)
% plot3(Z_coord{i},ksi_X{i},ksi_Y{i})
% end

