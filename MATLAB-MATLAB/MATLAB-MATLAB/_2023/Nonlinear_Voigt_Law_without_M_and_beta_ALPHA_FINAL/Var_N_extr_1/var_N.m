clc
clear
close all

%% Input data
m = 6;
h = 1/m;
% N = 30;
% Rod Parameters
d = 20*10^-3;           %[m]
l = 0.7;                %[m]
ro = 7800;              %[kg/m^3]
Elastic = 2e11;         %[Pa]

zeta_e = 0.025;         %[-]
zeta_V = 0.005;         %[-]
zeta_VV = 0.005;
alpha = 0.1;
mu_R = 1;
%%
[G00Int,G02Int,~,~] = get_matrix(m);
I=eye(m+1);        
E = eye(2);         
Z=zeros(m+1);       
R = [0 -1; 1 0]; 

%% Exclude the first and last nodes
[G00Int,G02Int,I,Z] = matrix_edging(G00Int,G02Int,I,Z,m);

%% Get matix M
M = get_matr_coeff_M(m);

%% Writing to the structure
[data] = Filling_structure(I,E,R,zeta_V,zeta_VV,zeta_e,G00Int,G02Int,mu_R,alpha,M,m);


%% Find crititcal velocity
N_start = 0;
step_N = 0.5;
Ncritical = get_N_crit (N_start,step_N,data);



%% %%% Solver Settings
T=[0,300];                       % the interval of the full study
opt=odeset('AbsTol',1e-6,'RelTol',1e-6);
x0 = zeros(2*(2*(m-1)),1); % vector of initial conditions
% Node number from 1 to (m-1)
number_node = 3;
x0(number_node+2) = 0.00005;
x0(number_node+3) = 0.00005;

%% %%% Get solve
dN = -0.1;
N_v = 24.5:dN:Ncritical+0.1;
end_time_ind = 20;

N_vector = [];
Amplitude = [];

for i=1:length(N_v)
    N_v(i)
   data.A0 = kron(I,(E-2*zeta_V*N_v(i)*R));
   data.A1 = 2*zeta_V*kron(I,E)+2*zeta_e*kron(G00Int,E);
   data.A2 = mu_R*kron(G00Int,E);
   data.N = N_v(i);
   [time,X] = SOLVE(data,T,x0,opt);
   [ksi_X_node] = get_KSI_X(X,number_node,m);
   N_vector(i) = N_v(i);
   Amplitude(i) = max(ksi_X_node);
   index = find(time>=time(end)-end_time_ind);
   [~,Extr_x] = ext(time(index),ksi_X_node(index));
   Extremum{i} = Extr_x;
   
   clear t X ksi_X_node Extr_x
end


%% Visualization N - Amplitude
figure
box on; grid on; hold on;
plot(N_vector,Amplitude,'.r','MarkerSize',15)


%% Visualization N - Extremum

figure
box on; grid on; hold on;
for i=1:length(Extremum)
    for j=1:length(Extremum{i})
        plot(N_vector(i),Extremum{i}(j),'.r','MarkerSize',15)
    end
end



function [t,X] = SOLVE(d,T,x0,opt)
    [t,X] = ode23t(@(t,X) solver(t,X,d.A0,d.A1,d.A2,d.G02Int,d.E,d.zeta_VV,d.alpha,d.M,d.N,d.R,d.m),T,x0,opt);
end

function [out] = get_KSI_X(X,number_node,m)
    for i=1:2:2*(m-1)
        ksi_X{i} = X(:,i);
    end
    ksi_X(:,2:2:end) = [];
    ksi_X = [zeros(length(ksi_X{1}),1),ksi_X,zeros(length(ksi_X{1}),1)];
    out = ksi_X{number_node+1};
end






