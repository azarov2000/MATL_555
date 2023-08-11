close all
clc
clear
%% Data
Gamma_0_vector = linspace(0,50,20);
alpha_vector = linspace(0,0.5,20);
zeta_VV_vector = linspace(0,0.05,40);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% Find Gamma_0  %%%%
for i=1:length(Gamma_0_vector)
    i
    Ampl_for_Gamma_0(i) = VAR(Gamma_0_vector(i),0,0,0.8);
end
%%
Gamma_0_RES = interp1(Ampl_for_Gamma_0,Gamma_0_vector,0.1);
%% Plot
figure;
box on; grid on; hold on;
plot(Gamma_0_vector,Ampl_for_Gamma_0,'.r','MarkerSize',16)
plot(Gamma_0_RES,0.1,'*b','MarkerSize',20)
ff = gca;
ff.FontSize = 16;
xlabel('\Gamma_{0}')
ylabel('Amplitude')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% Find alpha %%%%
for i=1:length(alpha_vector)
    i
    Ampl_for_alpha(i) = VAR(Gamma_0_RES,alpha_vector(i),0,0.8);
end 
%%
alpha_RES = interp1(Ampl_for_alpha,alpha_vector,0.095);
%% Plot
figure;
box on; grid on; hold on;
plot(alpha_vector,Ampl_for_alpha,'.r','MarkerSize',16)
plot(alpha_RES,0.095,'*b','MarkerSize',20)
ff = gca;
ff.FontSize = 16;
xlabel('\alpha')
ylabel('Amplitude')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% Find zeta_VV %%%%
for i=1:length(zeta_VV_vector)
    i
    Ampl_for_zeta_VV(i) = VAR(Gamma_0_RES,0,zeta_VV_vector(i),0.8);
end 
%%
zeta_VV_RES = interp1(Ampl_for_zeta_VV,zeta_VV_vector,0.095);
%% Plot
figure;
box on; grid on; hold on;
plot(zeta_VV_vector,Ampl_for_zeta_VV,'.r','MarkerSize',16)
plot(zeta_VV_RES,0.095,'*b','MarkerSize',20)
ff = gca;
ff.FontSize = 16;
xlabel('\eta_{ VV}')
ylabel('Amplitude')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% AFC %%%%
Freq_F_coeff_vector = 0.1:0.05:2;
for i=1:length(Freq_F_coeff_vector)
    i
%     Ampl(i) = VAR(Gamma_0_RES,alpha_RES,zeta_VV_RES,Freq_F_coeff_vector(i));
    Ampl(i) = VAR(10,0,0,Freq_F_coeff_vector(i));
end 
%% Plot
figure;
box on; grid on; hold on;
h1=stem(Freq_F_coeff_vector,Ampl);
set(get(h1,'BaseLine'),'LineStyle','-');
set(h1,'MarkerFaceColor','red');
ff = gca;
ff.FontSize = 16;
xlabel('\nu_{\Gamma} / \nu_{1}')
ylabel('Amplitude')


%%
function [MAX_AMPL] = VAR(Gamma_0,alpha,zeta_VV,F_coeff)
    %%
    % Cubic Feucht's law
    % Scheme: sealing-sealing
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

    N_z = 0;           % Осевая сила
    M_z = 0;           % Крутящий момент


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

    A0 = kron(I,(E-2*zeta_V*N*R)) - N_z*kron(G02Int,E) + M_z*kron(G01Int,R)*M;

    A1 = 2*zeta_V*kron(I,E)+2*zeta_e*kron(G00Int,E);

    A2 = mu_R*kron(G00Int,E);

    A2 = (A2+A2')./2;

    %% Writing to the structure
    [data] = Filling_structure(I,E,R,zeta_V,zeta_e,G00Int,G00,G01Int,G02Int,M,mu_R,M_z,N_z);
    %% Search for eigenvalues
    % [EigVal] = Calc_EigenValues(A0,A1,A2);

    %% Find Natural Frequency
    Number_freq = 3;
    dat.nat_freq = get_nat_freq(data,Number_freq);

    %% Поиск критической скорости

    % N_start = 0;
    % step_N = 0.5;
    % 
    % Ncritical = get_N_crit (N_start,step_N,data);


    %% Solution of the differential equation
    T=[0,15];                   % the interval of the full study
    opt=odeset('AbsTol',1e-6,'RelTol',1e-6);
%     opt = odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'OutputFcn', @(t, y, flag) odeoutput(t, y, flag, T(1), T(end) - T(1)));


    x0 = zeros(2*length(A0),1);  % vector of initial conditions
    % Node number from 1 to (m-1)
    number_node = round((m-1)/2);

    tic;
    [t,X] = ode23t(@(t,X) solver(t,X,A0,A1,A2,G02Int,G00,E,zeta_VV,alpha,M,N,R,m,Gamma_0,dat.nat_freq,F_coeff),T,x0,opt);
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
%%
    % Spectral analysis
%     dat.tNS = 5;
%     dat.number_node = round((m-1)/2);
%     dat.Power = 1;
%     
%     dat.ksi_x = ksi_X{dat.number_node+1};
%     dat.ksi_y = ksi_Y{dat.number_node+1};
%     dat.t = t;
%     dat.N = N;
%     
%     Spectral_analysis(dat)

    %% Nodes in 3D
    Number_end = 5;
    ind = find(t>T(end)-Number_end);
    % figure;
    % box on; grid on; hold on;
    % 
    % for i=1:length(Z_coord)
    % plot3(Z_coord{i}(ind),ksi_X{i}(ind),ksi_Y{i}(ind))
    % ff = gca;
    % ff.FontSize = 20;
    % 
    % xlabel('\zeta')
    % ylabel('\it x')
    % zlabel('\it y')
    % end

    % ylim([min(ksi_X{number_node+1}(ind)) max(ksi_X{number_node+1}(ind))]);
    % zlim([min(ksi_X{number_node+1}(ind)) max(ksi_X{number_node+1}(ind))]);

    % disp(['Max_disp = ', num2str(max(ksi_X{number_node+1}(ind)))])

    MAX_AMPL = max(ksi_X{number_node+1}(ind));

end


