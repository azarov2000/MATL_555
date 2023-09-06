close all
clc
clear
%% Кубический закон Фойхта
%% Заделка - Заделка

%% Входные параметры (для СЧ)
m = 60;
% zeta_e = 0.02;
zeta_e = 0.2;
zeta_V = 0.0001;
N_z = 0;
M_z = 0;

% Matrix and data creator
mc = matrix_creator(m,zeta_e,zeta_V,N_z,M_z);

% Find Natural Frequency
Number_freq = 4;
mc.natural_freq = get_nat_freq(mc,Number_freq);

% Find N_critical
% Ncritical = get_N_crit (0,0.1,mc);

%% Влияние момента на критичекую скорость
Moment_vector = linspace(-2,2,100);

for i=1:length(Moment_vector)
    tmp_mc = mc;
    tmp_mc.M_z = Moment_vector(i);
    N_crit_Mz(i) = get_N_crit (0,0.05,tmp_mc);
end
%%
figure;
box on; grid on; hold on;
plot(Moment_vector,N_crit_Mz,'.r','MarkerSize',14)

%% 2D - диаграмма (N and Moment)
N_vector = linspace(-10,10,100);
Moment_vector = linspace(-10,10,100);
figure;
box on; hold on;
for i=1:length(N_vector)
    i
   for j=1:length(Moment_vector)
       tmp_mc = mc;
       tmp_mc.N = N_vector(i);
       tmp_mc.M_z = Moment_vector(j);
       [A0,A1,A2] = Filling_matrix_A(tmp_mc);
       eig_values = polyeig(A0,A1,A2);
       if max(real(eig_values))>0
           markerColor = '.r';
           markerSize = 8;
       else
           markerColor = '.k';
           markerSize = 14;
       end 
       plot(N_vector(i),Moment_vector(j),markerColor,'MarkerSize',markerSize);
   end
end
xline(0,'--b')
yline(0,'--b')
xlabel('\Omega')
ylabel('M_{torsion}')
ff = gca;
ff.FontSize = 14;





%%
Moment_vector = linspace(0,1,1000);
tmp_mc = mc;
tmp_mc.N_z = 0;
tmp_mc.N = 2;
for i=1:length(Moment_vector)
    tmp_mc.M_z = Moment_vector(i);
    i
    [A0,A1,A2] = Filling_matrix_A(tmp_mc);
    DET_A2(i) = det(A2);
    
    eig_values = polyeig(A0,A1,A2);
    index = find(real(eig_values) == max(real(eig_values)));
    Max_RE_with_M(i) = real(eig_values(index(1)));
    Max_IM_with_M(i) = imag(eig_values(index(1)));
    
end
%%
figure
box on; grid on; hold on;
plot(Max_RE_with_M,Moment_vector,'-.r')
ff = gca;
ff.FontSize = 25;
xlabel("Max(Re(\lambda))")
ylabel("M_{ z}")

figure
box on; grid on; hold on;
plot(Max_RE_with_M,Max_IM_with_M,'-.r')
ff = gca;
ff.FontSize = 25;
xlabel("Max(Re(\lambda))")
ylabel("Max(Im(\lambda))")
%%
figure;
box on; grid on; hold on;
plot(Moment_vector,DET_A2,'-b')

%%
x0 = linspace(-2,2,50);
func = @(x) interp1(Moment_vector,DET_A2,x);

for i=1:length(x0)
    A(i) = fzero(func,x0(i));
end



