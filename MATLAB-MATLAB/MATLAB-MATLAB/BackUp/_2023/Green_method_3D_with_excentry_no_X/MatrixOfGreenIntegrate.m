    function [G00,G02,G03,G30]=MatrixOfGreenIntegrate(j,k,h)
    % G(r,s)
    
    r=h*(j-1); 
    s=h*(k-1);
    H=(r>s);
    %% Заделка - заделка
    C10=-1+3*s^2-2*s^3;     C12=6-12*s; C13=-12;
    C20=s-2*s^2+s^3;        C22=-4+6*s; C23=6; 
    C30=0;                  C32=0;      C33=0;
    C40=0;                  C42=0;      C43=0;

    G00 = H*(r-s)^3/6   +C10*r^3/6 +C20*r^2/2 +C30*r +C40;
    G02 = H*(r-s)       +C12*r^3/6 +C22*r^2/2 +C32*r +C42;
    G03 = -H            +C13*r^3/6 +C23*r^2/2 +C33*r +C43;
    G30 = H + C10;
    
    flag=h;
    if s==0 || s==1
        flag=flag/2;
    end
    G00=G00*flag;G02=G02*flag;G03=G03*flag;G30=G30*flag;
    end