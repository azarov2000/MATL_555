function [G00,G01,G10,G11]=MatrixOfGreen(j,k,h)

% G(r,s)
r=h*(j-1); 
s=h*(k-1);
H=(r>s);

        C10=-1+3*s^2-2*s^3;          C11=6*s-6*s^2;
        C20=s-2*s^2+s^3;             C21=1-4*s+3*s^2;
        C30=0;                       C31=0;
        C40=0;                       C41=0;

G00= H*(r-s)^3/6 +C10*r^3/6 +C20*r^2/2 +C30*r +C40;
G01=-H*(r-s)^2/2 +C11*r^3/6 +C21*r^2/2 +C31*r +C41;
G10= H*(r-s)^2/2 +C10*r^2/2 +C20*r     +C30;
G11=-H*(r-s)     +C11*r^2/2 +C21*r     +C31;

flag=1; % Нет домножения на 1/2 первого и последнего столбцов
G00=G00*flag;G01=G01*flag;G10=G10*flag;G11=G11*flag;