function [A0,A1,A2] = Filling_matrix_A(mc)

A0 = kron(mc.I,(mc.E-2*mc.zeta_V*mc.N*mc.R)) - ...
    mc.N_z*kron(mc.G02Int,mc.E) + mc.M_z*kron(mc.G01Int,mc.R)*mc.M;

A1 = 2*mc.zeta_V*kron(mc.I,mc.E)+2*mc.zeta_e*kron(mc.G00Int,mc.E);

A2 = kron(mc.G00Int,mc.E);

end

