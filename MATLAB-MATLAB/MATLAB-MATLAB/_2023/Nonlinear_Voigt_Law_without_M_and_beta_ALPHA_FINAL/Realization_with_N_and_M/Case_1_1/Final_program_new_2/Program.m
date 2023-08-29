close all
clc
clear
%% Кубический закон Фойхта
%% Заделка - Заделка

%% Входные параметры (для СЧ)
m = 12;
zeta_e = 0.02;
zeta_V = 0.0001;
N_z = 0;
M_z = 0;

% Matrix and data creator
mc = matrix_creator(m,zeta_e,zeta_V,N_z,M_z);

% Find Natural Frequency
Number_freq = 4;
mc.natural_freq = get_nat_freq(mc,Number_freq);

% Find N_critical
Ncritical = get_N_crit (0,0.1,mc);

%% Пямой ход
tmp_mc = mc;

N_vector_direct = 32:1:50;
T=[0,200];
opt=odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'OutputFcn', @(t, y, flag) odeoutput(t, y, flag, T(1), T(end) - T(1)));
tmp_mc.x0 = Initial_condition(tmp_mc,5e-3);
tmp_mc.F_coeff = 0;
tmp_mc.Number_end = 10;
tmp_mc.Base_number = 100;
tmp_mc.Gamma_0 = 0;
tmp_mc.alpha = 0.015;
tmp_mc.zeta_VV = 10^-5;
i = 1;
while i <= length(N_vector_direct)
    tmp_mc.N = N_vector_direct(i);
    [RES_N_dir{i}] = SOLVE(tmp_mc,T,opt);
    MAX_EXTR_DIR = RES_N_dir{i}.EXTR_1;
    Rel_tol_extr_dir = RES_N_dir{i}.Rel_tol_extr;
    i
    Rel_tol_extr_dir
    if Rel_tol_extr_dir <= 1e-3
        tmp_mc.x0 = RES_N_dir{i}.End_cond;
        T = [0, 200];
        i
        T_dir(i) = T(end);
        i = i+1;
    else 
        T(end) = T(end) + 200;
    end 
end

%% Ампилитуды
figure;
box on; grid on; hold on;
for i=1:length(RES_N_dir)
    plot(ones(length(RES_N_dir{i}.EXTR_1),1)*N_vector_direct(i),RES_N_dir{i}.EXTR_1,'.r','MarkerSize',12);
    ff = gca;
    ff.FontSize = 14;
    xlabel('N')
    ylabel('Amplitude')
end



%%
figure;
box on; grid on; hold on;
for i=1:length(RES_N_dir)
    plot(N_vector_direct(i),RES_N_dir{i}.MAX_AMPL_1,'.r','MarkerSize',15);
end
ff = gca;
ff.FontSize = 16;
xlabel('\Omega')
ylabel('Amplitude')
title(['m = ',num2str(m),'; \eta_{e} = ',num2str(zeta_e),...
       '; \eta_{V} = ',num2str(zeta_V),'; \Gamma_{0} = ',num2str(tmp_mc.Gamma_0),...
       '; \alpha = ',num2str(tmp_mc.alpha),'; \eta_{VV} = ',num2str(tmp_mc.zeta_VV)])

%% Обратный ход
tmp_mc = mc;

N_vector_reverce = 32:-0.2:10;
T = [0,800];
opt=odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'OutputFcn', @(t, y, flag) odeoutput(t, y, flag, T(1), T(end) - T(1)));
tmp_mc.x0 = Initial_condition(tmp_mc,5e-3);
tmp_mc.F_coeff = 0;
tmp_mc.Number_end = 2;
tmp_mc.Base_number = 50;
tmp_mc.Gamma_0 = 0;
tmp_mc.alpha = 0.015;
tmp_mc.zeta_VV = 10^-5;
i = 1;
while i <= length(N_vector_reverce)
    tmp_mc.N = N_vector_reverce(i);
    [RES_N_rev{i}] = SOLVE(tmp_mc,T,opt);
    MAX_EXTR_REV = RES_N_rev{i}.EXTR_1;
    Rel_tol_extr_rev = RES_N_rev{i}.Rel_tol_extr;
    if Rel_tol_extr_rev <= 1e-3
        tmp_mc.x0 = RES_N_rev{i}.End_cond;
        T = [0, 100];
        i
        T_rev(i) = T(end);
        i = i+1;
    else 
        T(end) = T(end) + 100;
    end 
end
%%
figure;
box on; grid on; hold on;
for i=1:33
    plot(N_vector_reverce(i),RES_N_rev{i}.MAX_AMPL_1,'.r','MarkerSize',15);
end
ff = gca;
ff.FontSize = 16;
xlabel('\Omega')
ylabel('Amplitude')
title(['m = ',num2str(m),'; \eta_{e} = ',num2str(zeta_e),...
       '; \eta_{V} = ',num2str(zeta_V),'; \Gamma_{0} = ',num2str(tmp_mc.Gamma_0),...
       '; \alpha = ',num2str(tmp_mc.alpha),'; \eta_{VV} = ',num2str(tmp_mc.zeta_VV)])