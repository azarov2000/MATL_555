function [A0,A1,A2]=test_matrix(N)
    
% Parameters
              %[-]
zeta_e = 0.025;     %[-]
zeta_V = 0.005;     %[-]


rootS = solver_for_alpha(1,1,20);
alpha1 = rootS(3);
E = eye(2); S=[0 1;-1 0];
A0 = alpha1^4 * (E+2*zeta_V*N*S);
A1 = 2*(alpha1^4*zeta_V+zeta_e)*E;
A2 = E;

end