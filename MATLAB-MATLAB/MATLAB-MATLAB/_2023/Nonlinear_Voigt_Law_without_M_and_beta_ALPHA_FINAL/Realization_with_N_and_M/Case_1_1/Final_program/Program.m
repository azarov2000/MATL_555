close all
clc
clear
%% Кубический закон Фойхта
%% Заделка - Заделка

%% Входные параметры (для СЧ)
m = 6;
N = 35;
zeta_e = 0.02;
zeta_V = 0.0001;
N_z = 0;
M_z = 0;

%% Matrix and data creator
mc = matrix_creator(m,N,zeta_e,zeta_V,N_z,M_z);

 %% Find Natural Frequency
Number_freq = 4;
mc.natural_freq = get_nat_freq(mc,Number_freq);





