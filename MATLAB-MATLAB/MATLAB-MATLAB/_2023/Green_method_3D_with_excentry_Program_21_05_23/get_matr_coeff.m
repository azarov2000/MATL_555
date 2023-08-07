function [dat] = get_matr_coeff(dat)

a_psi_0 = kron(dat.I,dat.E+2*dat.zeta_V*dat.N*dat.S)-dat.Mom*kron(dat.G03Int,dat.S);
a_X_0 = -12*dat.Mom*kron(dat.Fg,dat.E);

b_psi_0 = 12*kron(dat.I(end,:),dat.E);
b_X_0 = -12*dat.Mom*dat.Fg_ticks*dat.S+12*(dat.E+2*dat.zeta_V*dat.N*dat.S)+6*dat.Mom*dat.S+dat.kappa_B*dat.E;

A0 = [a_psi_0, a_X_0;
      b_psi_0, b_X_0];


a_psi_1 = 2*dat.zeta_V*kron(dat.I,dat.E)-2*dat.zeta_e*kron(dat.G00Int,dat.E)-2*dat.N*dat.beta*kron(dat.G02Int,dat.S);
a_X_1 = -2*dat.zeta_e*kron(dat.F0,dat.E)-2*dat.beta*dat.N*kron(dat.F2,dat.S);

b_psi_1 = -2*dat.zeta_e*kron(dat.G30Int(end,:),dat.E)-2*dat.N*dat.beta*kron(dat.T',dat.S);
b_X_1 = -2*dat.zeta_e*dat.F0_ticks*dat.E-2*dat.beta*dat.N*dat.F2_ticks*dat.S+24*dat.zeta_V*dat.E+2*dat.eta_B*dat.E;
  
A1 = [a_psi_1, a_X_1;
      b_psi_1, b_X_1];
  
a_psi_2 = kron(dat.G00Int,dat.E)-dat.beta*kron(dat.G02Int,dat.E);
a_X_2 = kron(dat.F0,dat.E)-dat.beta*kron(dat.F2,dat.E);
b_psi_2 = kron(dat.G30Int(end,:),dat.E)-dat.beta*kron(dat.T',dat.E);
b_X_2 = dat.F0_ticks*dat.E-dat.beta*dat.F2_ticks*dat.E+dat.mu_B*dat.E; 
 
A2 = [a_psi_2, a_X_2;
      b_psi_2, b_X_2];
 
g_psi = kron(dat.G00Int,dat.E);
g_X = kron(dat.G30Int(end,:),dat.E);

G = [g_psi, zeros(length(kron(dat.G00Int,dat.E)));
    zeros(2,length(kron(dat.G00Int,dat.E))) ,g_X];


dat.A0 = A0;
dat.A1 = A1;
dat.A2 = A2;
dat.G = G;

end

