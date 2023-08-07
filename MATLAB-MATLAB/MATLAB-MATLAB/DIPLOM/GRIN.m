function [G00,G0s,Gr0,Grs]=GRIN(r,s,nameS)

H=(r>s);
switch nameS
    case 1      %    C30=C40=0-------Clamped-Hinged
        C10=(-2+3*s^2-s^3)/2;        C1s=(6*s-3*s^2)/2;
        C20=(2*s-3*s^2+s^3)/2;       C2s=(2-6*s+3*s^2)/2;
        C30=0;                       C3s=0;
        C40=0;                       C4s=0;

    case 2      %    C20=C40=0-------Hinged-Hinged
        C10=s-1;                     C1s=1;
        C20=0;                       C2s=0;
        C30=(2*s-3*s^2+s^3)/6;       C3s=(2-6*s+3*s^2)/6;
        C40=0;                       C4s=0;
        
    case 3      %    C30=C40=0--------Clamped-Free
        C10=-1;                      C1s=0;
        C20=s;                       C2s=1;
        C30=0;                       C3s=0;
        C40=0;                       C4s=0;

    case 4      %    C30=C40=0--------Clamped-Clamped
        C10=-1+3*s^2-2*s^3;          C1s=6*s-6*s^2;
        C20=s-2*s^2+s^3;             C2s=1-4*s+3*s^2;
        C30=0;                       C3s=0;
        C40=0;                       C4s=0;
               
end

G00= H*(r-s)^3/6 +C10*r^3/6 +C20*r^2/2 +C30*r +C40;
G0s=-H*(r-s)^2/2 +C1s*r^3/6 +C2s*r^2/2 +C3s*r +C4s;
Gr0= H*(r-s)^2/2 +C10*r^2/2 +C20*r     +C30;
Grs=-H*(r-s)     +C1s*r^2/2 +C2s*r     +C3s;