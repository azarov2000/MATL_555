close all
clc
clear
%% 
n = 1;
% Parameters
d = 20*10^-3;       %[m]
l = 0.7;            %[m]
ro = 7800;          %[kg/m^3]

N = 22;             %[-]
eps_d = d/l;        %[-]
beta = (eps_d^2)/16;%[-]
zeta_e = 0.025;     %[-]
zeta_V = 0.005;     %[-]


%% Argan diagram and critical_a
N_start = 0;
N_step = 1;
v_zeta_V = linspace(0.002,0.1,30);
for i=1:length(v_zeta_V)
    N_cr(i) = Critical_Speed_Finder_galerkin(N_start,N_step,n,beta,v_zeta_V(i),zeta_e);
end 

%% Variation of internal damping

figure;
box on; grid on; hold on; 
plot(v_zeta_V,N_cr_galerkin,'.-','MarkerSize',12)




















