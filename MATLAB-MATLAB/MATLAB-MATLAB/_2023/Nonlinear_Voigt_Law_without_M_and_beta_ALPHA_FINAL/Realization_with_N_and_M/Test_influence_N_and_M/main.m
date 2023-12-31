close all
clc
clear
% Cubic Feucht's law
% Scheme: sealing-sealing
warning('on')
m = 12;
h = 1/m;
N = 25;                  % 0 

% Rod Parameters
d = 20*10^-3;           %[m]
l = 0.7;                %[m]
ro = 7800;              %[kg/m^3]
Elastic = 2e11;         %[Pa]

zeta_e = 0.025;         %[-]    - 0
zeta_V = 0.005;         %[-]    - 0 


zeta_VV = 0.005;
alpha = 0.1;
mu_R = 1;

N_z = 0;             % Осевая сила     0
M_z = 0.1;           % Крутящий момент  0

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

%% Writing to the structure
[data] = Filling_structure(I,E,R,zeta_V,zeta_e,G00Int,G01Int,G02Int,M,mu_R,M_z,N_z);
%% Search for eigenvalues

[EigVal] = Calc_EigenValues(A0,A1,A2);

%% Find_critical_velocity
N_start = 0;
step_N = 0.5;

Ncritical = get_N_crit (N_start,step_N,data);


%% Influence tors_Moment on critical velocirty
M_z_vector = linspace(-10,10,10000);
for i=1:length(M_z_vector)
    i
    data_tmp_M = data;
    data_tmp_M.M_z = M_z_vector(i);
    N_crit_M(i) = get_N_crit(0.01,0.1,data_tmp_M);
end
%%
figure
box on; hold on; grid on;
plot(M_z_vector, N_crit_M,'.r','MarkerSize',15)
ff = gca;
ff.FontSize = 25;
xlabel("M_{ z}")
ylabel("N_{ crit}")
%% Influence axis_Force on critical velocirty
% N_z_vector = linspace(-10,10,1000);
% for i=1:length(N_z_vector)
%     i
%     data_tmp_N = data;
%     data_tmp_N.N_z = N_z_vector(i);
%     N_crit_N(i) = get_N_crit(0.01,0.1,data_tmp_N);
% end
% %%
% figure
% box on; hold on; grid on;
% plot(N_z_vector, N_crit_N,'.r','MarkerSize',15)
% xlabel("N_{ z}")
% ylabel("N_{ crit}")


%% max(Re(lyambda)) x M_z
M_z_vector = linspace(-10,10,10000);
dd = data;
dd.N_z = 0;
dd.N = 0;
for i=1:length(M_z_vector)
    dd.M_z = M_z_vector(i);
    
    dd.A0 = kron(I,(E-2*zeta_V*dd.N*R)) - dd.N_z*kron(G02Int,E) + dd.M_z*kron(G01Int,R)*M;
    dd.A1 = 2*zeta_V*kron(I,E)+2*zeta_e*kron(G00Int,E);
    dd.A2 = mu_R*kron(G00Int,E);
    
    eig_values = polyeig(dd.A0,dd.A1,dd.A2);
    index = find(real(eig_values) == max(real(eig_values)));
    Max_RE_with_M(i) = real(eig_values(index(1)));
    Max_IM_with_M(i) = imag(eig_values(index(1)));
    
end

figure
box on; grid on; hold on;
plot(Max_RE_with_M,M_z_vector)
ff = gca;
ff.FontSize = 25;
xlabel("Max(Re(\lambda))")
ylabel("M_{ z}")

figure
box on; grid on; hold on;
plot(Max_RE_with_M,Max_IM_with_M)
ff = gca;
ff.FontSize = 25;
xlabel("Max(Re(\lambda))")
ylabel("Max(Im(\lambda))")

%% max(Re(lyambda)) x N_z

N_z_vector = linspace(-100,100,10000);
dd = data;
dd.M_z = 0;
dd.N = 0;
for i=1:length(N_z_vector)
    dd.N_z = N_z_vector(i);
    
    dd.A0 = kron(I,(E-2*zeta_V*dd.N*R)) - dd.N_z*kron(G02Int,E) + dd.M_z*kron(G01Int,R)*M;
    dd.A1 = 2*zeta_V*kron(I,E)+2*zeta_e*kron(G00Int,E);
    dd.A2 = mu_R*kron(G00Int,E);
    
    eig_values = polyeig(dd.A0,dd.A1,dd.A2);
    index = find(real(eig_values) == max(real(eig_values)));
    Max_RE_with_M(i) = real(eig_values(index(1)));
    Max_IM_with_M(i) = imag(eig_values(index(1)));
    
end

figure
box on; grid on; hold on;
plot(Max_RE_with_M,N_z_vector)
ff = gca;
ff.FontSize = 25;
xlabel("Max(Re(\lambda))")
ylabel("N_{ z}")

figure
box on; grid on; hold on;
plot(Max_RE_with_M,Max_IM_with_M)
ff = gca;
ff.FontSize = 25;
xlabel("Max(Re(\lambda))")
ylabel("Max(Im(\lambda))")
%% Solution of the differential equation

opt=odeset('AbsTol',1e-6,'RelTol',1e-6);

T=[0,600];                   % the interval of the full study
x0 = zeros(2*length(A0),1);  % vector of initial conditions
% Node number from 1 to (m-1)
number_node = 3;

x0(number_node+2) = 0.00005;
x0(number_node+3) = 0.00005;


tic;
[t,X] = ode23t(@(t,X) solver(t,X,A0,A1,A2,G02Int,E,zeta_VV,alpha,M,N,R,m),T,x0,opt);
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
figure;
box on; grid on; hold on;

for i=1:length(Z_coord)
plot3(Z_coord{i},ksi_X{i},ksi_Y{i})
end

