function [G00Int,G01Int,G02Int,G00,G01,G02,I,Z] = matrix_edging(G00Int,G01Int,G02Int,G00,G01,G02,I,Z,m)

Lim = m;
G00Int = G00Int(2:Lim,2:Lim);
G01Int = G01Int(2:Lim,2:Lim);
G02Int = G02Int(2:Lim,2:Lim);

G00 = G00(2:Lim,2:Lim);
G01 = G01(2:Lim,2:Lim);
G02 = G02(2:Lim,2:Lim);

I = I(2:Lim,2:Lim);
Z = Z(2:Lim,2:Lim);

end




