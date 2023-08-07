function [G00Int,G01Int,G02Int,G10Int,G11Int,G12Int,...
          G00,G01,G02,G10,G11,G12,I,EC,EC_1,Z]=matrix_edging(G00Int,G01Int,G02Int,G10Int,G11Int,G12Int,...
          G00,G01,G02,G10,G11,G12,I,EC,EC_1,Z)


G00Int = G00Int(2:end,2:end);
G01Int = G01Int(2:end,2:end);
G02Int = G02Int(2:end,2:end);
G10Int = G10Int(2:end,2:end);
G11Int = G11Int(2:end,2:end);
G12Int = G12Int(2:end,2:end);

G00 = G00(2:end,2:end);
G01 = G01(2:end,2:end);
G02 = G02(2:end,2:end);
G10 = G10(2:end,2:end);
G11 = G11(2:end,2:end);
G12 = G12(2:end,2:end);

I = I(2:end,2:end);
EC = EC(2:end,2:end);
EC_1 = EC_1(2:end,2:end);
Z = Z(2:end,2:end);


end




