function [G00,G01,G10,G11,G00Int,G01Int,G10Int,G11Int,I,Z,EC,m,jC,DE,mnew] = Matrix_Edging(G00,G01,G10,G11,G00Int,G01Int,G10Int,G11Int,I,Z,EC,DE,m,jC,L)
    
         G00=G00(2:L-1,2:L-1); G01=G01(2:L-1,2:L-1); G10=G10(2:L-1,2:L-1);
         G11=G11(2:L-1,2:L-1); G00Int=G00Int(2:L-1,2:L-1); G01Int=G01Int(2:L-1,2:L-1);
         G10Int=G10Int(2:L-1,2:L-1); G11Int=G11Int(2:L-1,2:L-1);
         I=I(2:L-1,2:L-1); Z=Z(2:L-1,2:L-1); EC=EC(2:L-1,2:L-1);
         DE=DE(2:L-1,2:L-1);
         mnew=m-2;  % Стало на 2 элемента меньше
         jC=jC-1;   % Смешение индексации (т.к. убираем первый узел)
end