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

%% Поиск значения Г_0 (приложение силы к срединному узлу)
tmp_mc = mc;

Gamma_0_vector = linspace(0,10,2);
T=[0,200];
opt=odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'OutputFcn', @(t, y, flag) odeoutput(t, y, flag, T(1), T(end) - T(1)));
tmp_mc.x0 = Initial_condition(tmp_mc,0);
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
tmp_mc.x0 = Initial_condition(tmp_mc,0);
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

zeta_VV_vector = [0:0.00005:0.005,0.01:0.005:0.1];
T=[0,100];
opt=odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'OutputFcn', @(t, y, flag) odeoutput(t, y, flag, T(1), T(end) - T(1)));
tmp_mc.x0 = Initial_condition(tmp_mc,0);
tmp_mc.N = 0;
tmp_mc.F_coeff = 0.8;
tmp_mc.Number_end = 5;
tmp_mc.Gamma_0 = Gamma_0_RES;
% tmp_mc.alpha = 0.015;
tmp_mc.alpha = 0;

for i=1:length(zeta_VV_vector)
    tmp_mc.zeta_VV = zeta_VV_vector(i);
    [t{i},ksi_X{i},ksi_Y{i},MAX_AMPL_zeta_VV(i)] = SOLVE(tmp_mc,T,opt);
end
%%
zeta_VV_RES = interp1(MAX_AMPL_zeta_VV,zeta_VV_vector,0.095);
tmp_zeta_VV = interp1(zeta_VV_vector,MAX_AMPL_zeta_VV,1e-5);
%%
figure;
box on; grid on; hold on;
plot(zeta_VV_vector,MAX_AMPL_zeta_VV,'.r','MarkerSize',20)
plot(zeta_VV_RES,0.095,'*k','MarkerSize',20)
plot(1e-5,tmp_zeta_VV,'*b','MarkerSize',20)
yline(0.095+0.095*0.10,'--r','В пределах 10%','FontSize',14,'LineWidth',2)
yline(0.095-0.095*0.10,'--r','LineWidth',2)

ff = gca;
ff.FontSize = 20;
title(['m = ',num2str(m),'; \eta_{e} = ',num2str(zeta_e),...
       '; \eta_{V} = ',num2str(zeta_V),'; \Gamma_{0} = ',num2str(tmp_mc.Gamma_0),...
       '; \alpha = ',num2str(tmp_mc.alpha),'; \Omega = ',num2str(tmp_mc.N)])
xlabel('\eta_{VV}');
ylabel('Amplitude');
% legend("AutoUpdate","on")

%% Реализации
zeta_VV_number = 20;
figure;
box on; grid on; hold on;
plot(t{zeta_VV_number},ksi_X{zeta_VV_number}{7})
ff = gca;
ff.FontSize = 20;
title(['m = ',num2str(m),'; \eta_{e} = ',num2str(zeta_e),...
       '; \eta_{V} = ',num2str(zeta_V),'; \alpha = ',num2str(tmp_mc.alpha),...
       '; \eta_{VV} = ',num2str(zeta_VV_vector(zeta_VV_number)),'; \Omega = ',num2str(tmp_mc.N)])
xlabel('t');
ylabel('\xi_{x}');


%% Построрение АЧХ (приложение силы к срединному узлу)
tmp_mc = mc;

F_coeff_vector = 0.5:0.02:1.5;
T=[0,100];
opt=odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'OutputFcn', @(t, y, flag) odeoutput(t, y, flag, T(1), T(end) - T(1)));
tmp_mc.x0 = Initial_condition(tmp_mc,0);
tmp_mc.N = 0;
tmp_mc.Number_end = 5;
tmp_mc.Gamma_0 = Gamma_0_RES;
tmp_mc.alpha = 0.015;
tmp_mc.zeta_VV = 10^-5;

for i=1:length(F_coeff_vector)
    tmp_mc.F_coeff = F_coeff_vector(i);
    [t{i},ksi_X{i},ksi_Y{i},MAX_AMPL(i)] = SOLVE(tmp_mc,T,opt);
end
%%
figure;
box on; grid on; hold on;
h1=stem(F_coeff_vector,MAX_AMPL);
set(get(h1,'BaseLine'),'LineStyle','-');
set(h1,'MarkerFaceColor','red');
ff = gca;
ff.FontSize = 16;
xlabel('\nu_{\Gamma} / \nu_{1}')
ylabel('Amplitude')
title(['m = ',num2str(m),'; \eta_{e} = ',num2str(zeta_e),...
       '; \eta_{V} = ',num2str(zeta_V),'; \Gamma_{0} = ',num2str(tmp_mc.Gamma_0),...
       '; \alpha = ',num2str(tmp_mc.alpha),'; \eta_{VV} = ',num2str(tmp_mc.zeta_VV),...
       '; \Omega = ',num2str(tmp_mc.N)])


%% Пямой ход и обратный
tmp_mc = mc;

N_vector_direct = 33:1:50;
T=[0,2000];
opt=odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'OutputFcn', @(t, y, flag) odeoutput(t, y, flag, T(1), T(end) - T(1)));
tmp_mc.x0 = Initial_condition(tmp_mc,5e-3);
tmp_mc.F_coeff = 0;
tmp_mc.Number_end = 5;
tmp_mc.Gamma_0 = 0;
tmp_mc.alpha = 0.015;
tmp_mc.zeta_VV = 10^-5;
i = 1;
while i <= length(N_vector_direct)
    tmp_mc.N = N_vector_direct(i);
    [RES_N_dir{i}] = SOLVE(tmp_mc,T,opt);
    disp(RES_N_dir{i}.Rel_tol_extr)
    if RES_N_dir{i}.Rel_tol_extr <= 1e-8
        tmp_mc.x0 = RES_N_dir{i}.End_cond;
        T = [0, 200];
        i
        T_res(i) = T(end);
        i = i+1;
    else 
        T(end) = T(end) + 200;
    end 
end


% for i=1:length(N_vector_direct)
%     tmp_mc.N = N_vector_direct(i);
%     [RES_N_dir{i}] = SOLVE(tmp_mc,T,opt);
%     if i==1
%         T(end) = T(end);
%     end
% end

%%
figure;
box on; grid on; hold on;
h1=stem(N_vector_direct,MAX_AMPL_direct_N);
set(get(h1,'BaseLine'),'LineStyle','-');
set(h1,'MarkerFaceColor','red');
ff = gca;
ff.FontSize = 16;
xlabel('\Omega')
ylabel('Amplitude')
title(['m = ',num2str(m),'; \eta_{e} = ',num2str(zeta_e),...
       '; \eta_{V} = ',num2str(zeta_V),'; \Gamma_{0} = ',num2str(tmp_mc.Gamma_0),...
       '; \alpha = ',num2str(tmp_mc.alpha),'; \eta_{VV} = ',num2str(tmp_mc.zeta_VV)])

