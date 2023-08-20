close all
clc
clear
%% Кубический закон Фойхта
%% Заделка - Заделка

%% Входные параметры (для СЧ)
m = 12;
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

%% Find N_critical
Ncritical = get_N_crit (0,0.1,mc);

%% Поиск значения Г_0 (приложение силы к срединному узлу)
tmp_mc = mc;

Gamma_0_vector = linspace(0,10,2);
T=[0,200];
opt=odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'OutputFcn', @(t, y, flag) odeoutput(t, y, flag, T(1), T(end) - T(1)));
tmp_mc.N = 0;
tmp_mc.zeta_VV = 0;
tmp_mc.alpha = 0;
tmp_mc.F_coeff = 0.8;
tmp_mc.Number_end = 5;

for i=1:length(Gamma_0_vector)
    tmp_mc.Gamma_0 = Gamma_0_vector(i);
    [t{i},ksi_X{i},ksi_Y{i},MAX_AMPL_Gamma(i)] = SOLVE(tmp_mc,T,opt);
end
Gamma_0_RES = interp1(MAX_AMPL_Gamma,Gamma_0_vector,0.1);
%%
figure;
box on; grid on; hold on;
plot(Gamma_0_vector,MAX_AMPL_Gamma,'.r','MarkerSize',20)
plot(Gamma_0_RES,0.1,'*k','MarkerSize',20)
ff = gca;
ff.FontSize = 20;
title(['m = ',num2str(m),'; \eta_{e} = ',num2str(zeta_e),...
       '; \eta_{V} = ',num2str(zeta_V),'; \alpha = ',num2str(tmp_mc.alpha),...
       '; \eta_{VV} = ',num2str(tmp_mc.zeta_VV),'; \Omega = ',num2str(tmp_mc.N)])
xlabel('\Gamma_{0}');
ylabel('Amplitude');
%% Реализации
figure;
box on; grid on; hold on;
plot(t{2},ksi_X{2}{7})
ff = gca;
ff.FontSize = 20;
title(['m = ',num2str(m),'; \eta_{e} = ',num2str(zeta_e),...
       '; \eta_{V} = ',num2str(zeta_V),'; \alpha = ',num2str(tmp_mc.alpha),...
       '; \eta_{VV} = ',num2str(tmp_mc.zeta_VV),'; \Omega = ',num2str(tmp_mc.N)])
xlabel('t');
ylabel('\xi_{x}');

%% Поиск значения alpha (приложение силы к срединному узлу)
tmp_mc = mc;

alpha_vector = 0:0.005:0.1;
T=[0,100];
opt=odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'OutputFcn', @(t, y, flag) odeoutput(t, y, flag, T(1), T(end) - T(1)));
tmp_mc.N = 0;
tmp_mc.zeta_VV = 0;
tmp_mc.F_coeff = 0.8;
tmp_mc.Number_end = 5;
tmp_mc.Gamma_0 = Gamma_0_RES;

for i=1:length(alpha_vector)
    tmp_mc.alpha = alpha_vector(i);
    [t{i},ksi_X{i},ksi_Y{i},MAX_AMPL_alpha(i)] = SOLVE(tmp_mc,T,opt);
end
alpha_RES = interp1(MAX_AMPL_alpha,alpha_vector,0.095);
%%
figure;
box on; grid on; hold on;
plot(alpha_vector,MAX_AMPL_alpha,'.r','MarkerSize',20)
plot(alpha_RES,0.095,'*k','MarkerSize',20)
yline(0.095+0.095*0.01,'--r','В пределах 1%','FontSize',14,'LineWidth',2)
yline(0.095-0.095*0.01,'--r','LineWidth',2)

ff = gca;
ff.FontSize = 20;
title(['m = ',num2str(m),'; \eta_{e} = ',num2str(zeta_e),...
       '; \eta_{V} = ',num2str(zeta_V),'; \Gamma_{0} = ',num2str(tmp_mc.Gamma_0),...
       '; \eta_{VV} = ',num2str(tmp_mc.zeta_VV),'; \Omega = ',num2str(tmp_mc.N)])
xlabel('\alpha');
ylabel('Amplitude');

%% Реализации
alpha_number = 4;
figure;
box on; grid on; hold on;
plot(t{alpha_number},ksi_X{alpha_number}{7})
ff = gca;
ff.FontSize = 20;
title(['m = ',num2str(m),'; \eta_{e} = ',num2str(zeta_e),...
       '; \eta_{V} = ',num2str(zeta_V),'; \alpha = ',num2str(alpha_vector(alpha_number)),...
       '; \eta_{VV} = ',num2str(tmp_mc.zeta_VV),'; \Omega = ',num2str(tmp_mc.N)])
xlabel('t');
ylabel('\xi_{x}');


%% Поиск значения ета_VV (приложение силы к срединному узлу)
tmp_mc = mc;

zeta_VV_vector = 0:0.005:0.1;
T=[0,100];
opt=odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'OutputFcn', @(t, y, flag) odeoutput(t, y, flag, T(1), T(end) - T(1)));
tmp_mc.N = 0;
tmp_mc.F_coeff = 0.8;
tmp_mc.Number_end = 5;
tmp_mc.Gamma_0 = Gamma_0_RES;
tmp_mc.alpha = 0.015;

for i=1:length(zeta_VV_vector)
    tmp_mc.zeta_VV = zeta_VV_vector(i);
    [t{i},ksi_X{i},ksi_Y{i},MAX_AMPL_zeta_VV(i)] = SOLVE(tmp_mc,T,opt);
end
zeta_VV_RES = interp1(MAX_AMPL_zeta_VV,zeta_VV_vector,0.095);








