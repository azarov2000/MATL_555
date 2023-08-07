close all
clc
clear
%%
% Parameters
n = 2;
d = 20*10^-3;       %[m]
l = 0.7;            %[m]
ro = 7800;          %[kg/m^3]

N = 12;             %[-]
eps_d = d/l;        %[-]
beta = (eps_d^2)/16;%[-]
zeta_e = 0.025;     %[-]
zeta_V = 0.005;     %[-]

zeta_Dksi = 0.0025;  % [-] - коэффициент внешнего линейного демфирования диска 
zeta_Dte = 0.05;     % [-] - коэффициент внешнего углового демфирования диска

type_BC = 1;        %[-] тип

zeta_C = 0.5;

%%
roots_alpha = solver_for_alpha(type_BC,1,30);
%%

E = eye(2); I = eye(n); S=[0 1;-1 0];
alpha_vect = roots_alpha(1:n);
for i=1:n
    ALPHA(i,i) = ((alpha_vect(i)))^4;
end

% Obtaining orthonormal psi functions
inc_z = 1000;
z = linspace(0,1,inc_z);
dz = mean(diff(z));
for i=1:n   % index alpha
    K1{i} = @(z) krylovF(z,alpha_vect(i),1);
    K2{i} = @(z) krylovF(z,alpha_vect(i),2);
    K3{i} = @(z) krylovF(z,alpha_vect(i),3);
    K4{i} = @(z) krylovF(z,alpha_vect(i),4);
    Psi_tmp{i} = @(z) (K3{i}(z)-(K3{i}(1)/K4{i}(1))*K4{i}(z));  
    % Normalization
    Norm_C(i) = trapz(z,(Psi_tmp{i}(z)).^2);
    Psi{i} = @(z) Psi_tmp{i}(z)/sqrt(Norm_C(i));   
end
% Diff 1 and 2
for i=1:n
    D_test{i} = @(z) (Norm_C(i))^(-0.5)*(alpha_vect(i))*(K2{i}(z)-(K3{i}(1)/K4{i}(1))*K3{i}(z));
    DD_test{i} = @(z) (Norm_C(i))^(-0.5)*(alpha_vect(i))^2*(K1{i}(z)-(K3{i}(1)/K4{i}(1))*K2{i}(z));
end

%% Filling in the matrix G

for i=1:n
    for k=1:n
        G(i,k) = trapz(z,Psi{i}(z).*DD_test{k}(z));
        F(i,k) = Psi{i}(zeta_C).*Psi{k}(zeta_C);
        H(i,k) = D_test{i}(zeta_C).*D_test{k}(zeta_C);
    end
end

% All matrix
A0 = kron(ALPHA,E)+2*zeta_V*N*kron(ALPHA,S);
A1 = -2*beta*N*kron(G,S) + kron(2*(zeta_V*ALPHA+zeta_e*I),E) + 2*zeta_Dksi*kron(F,E)+...
    2*N*beta*kron(H,S)+2*zeta_Dte*kron(H,E);
A2 = -beta*kron(G,E)+kron(I,E)+kron(F,E)+beta*kron(H,E);


%% 
EigenValues = polyeig(A0,A1,A2)


