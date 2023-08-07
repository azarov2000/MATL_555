function [G00Int,G01Int,G02Int,G20Int,G21Int,...
          G00,G01,G02,G20,G21,J,I,EC,Z]=matrix_edging(G00Int,G01Int,G02Int,G20Int,G21Int,...
          G00,G01,G02,G20,G21,J,I,EC,Z,m)

Lim = m;
G00Int = G00Int(2:Lim,2:Lim);
G01Int = G01Int(2:Lim,2:Lim);
G02Int = G02Int(2:Lim,2:Lim);
G20Int = G20Int(2:Lim,2:Lim);
G21Int = G21Int(2:Lim,2:Lim);

G00 = G00(2:Lim,2:Lim);
G01 = G01(2:Lim,2:Lim);
G02 = G02(2:Lim,2:Lim);
G20 = G20(2:Lim,2:Lim);
G21 = G21(2:Lim,2:Lim);

J = J(2:Lim,2:Lim);
I = I(2:Lim,2:Lim);
EC = EC(2:Lim,2:Lim);
Z = Z(2:Lim,2:Lim);


end




