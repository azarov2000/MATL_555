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

%% Critical Speed
N_start = 0;
N_step = 1; 
N_cr_1 = Critical_Speed_Finder_galerkin(N_start,N_step,1,beta,zeta_V,zeta_e)
N_cr_2 = Critical_Speed_Finder_galerkin(N_start,N_step,2,beta,zeta_V,zeta_e)


%% Argan diagram and critical_a
N_start = 0;
N_step = 1;
v_zeta_V = linspace(0.002,0.1,50);
for i=1:length(v_zeta_V)
    N_cr(i) = Critical_Speed_Finder_galerkin(N_start,N_step,n,beta,v_zeta_V(i),zeta_e);
end 
%%
N_cr_analytics = @(zeta_V) pi^2*(1+zeta_e/(pi^4*zeta_V))*sqrt(1/(1+pi^2*beta-2*pi^2*(1+zeta_e/(pi^4*zeta_V))*beta));
%% Variation of internal damping
figure;
box on; grid on; hold on; 
plot(v_zeta_V,N_cr,'.-','MarkerSize',12)
for i=1:length(v_zeta_V)
    plot(v_zeta_V(i),N_cr_analytics(v_zeta_V(i)),'.r','MarkerSize',12)
end
%% Проверка

Omega1 = @(N) (pi^2*(N*beta+sqrt(1+beta*(pi^2+N^2*beta))))/(1+pi^2*beta);
Omega2 = @(N) (pi^2*(N*beta-sqrt(1-beta*(pi^2+N^2*beta))))/(1+pi^2*beta);
N_v = linspace(0,20,100);
figure;
box on; grid on; hold on;
for i=1:length(N_v)
    plot(N_v(i),Omega1(N_v(i)),'.r','MarkerSize',12)
    plot(N_v(i),Omega2(N_v(i)),'.r','MarkerSize',12)
end

%eig = @(N) polyeig(pi^4*(1-2*i*zeta_V*N),2*(pi^4*zeta_V+zeta_e-pi^2*beta*2*i*N),1+pi^2*beta)

eig = @(N) polyeig(pi^4*(1-2*i*zeta_V*N),2*(pi^4*zeta_V+zeta_e-pi^2*beta*2*i*N),1+pi^2*beta)



















