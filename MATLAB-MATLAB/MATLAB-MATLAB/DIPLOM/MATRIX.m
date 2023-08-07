function [A2,A1,A0]=MATRIX(z_C,N)

global beta zeta_e zeta_i Zeta_e Zeta_te nameS

[G00,G0s,Gr0,Grs]=GRIN(z_C,z_C,nameS);

Z=(zeta_e+zeta_i+Zeta_e);Zt=2*beta*Zeta_te;

O=zeros(2,2);E=eye(2,2);S=[0 1;-1 0];
A0=[E+2*zeta_i*N*G00*S, O;...                             
    2*zeta_i*N*Gr0*E,   E];

A1=[2*Z*G00*E,  2*G0s*(-beta*N*E+Zt*S);...
    -2*Z*Gr0*S, 2*Grs*(Zt*E+beta*N*S)];

A2=[G00*E,  G0s*beta*S;...
    -Gr0*S, beta*Grs*E];