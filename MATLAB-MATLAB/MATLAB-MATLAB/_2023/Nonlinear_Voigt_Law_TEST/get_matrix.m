function [G00Int,G01Int,G02Int,G03Int,G20Int,G21Int,...
          G00,G01,G02,G03,G20,G21,J]=get_matrix(m)
    
    h = 1/m;
    
% Вычисление матриц функций
    for j=1:m+1     % проходим по строкам матриц
        for k=1:m+1 % проходим по столбцам матриц
            [G00Int(j,k),G01Int(j,k),G02Int(j,k),G03Int(j,k),G20Int(j,k),G21Int(j,k)]=MatrixOfGreenIntegrate(j,k,h);
            [G00(j,k),G01(j,k),G02(j,k),G03(j,k),G20(j,k),G21(j,k)]=MatrixOfGreen(j,k,h);
            [J(j,k)]=MatrixJ(j,k,h);
        end
    end
    
    function [G00,G01,G02,G03,G20,G21]=MatrixOfGreenIntegrate(j,k,h)

    % G(r,s)
    r=h*(j-1); 
    s=h*(k-1);
    H=(r>=s);
    
    %% Заделка - заделка
    C10=-1+3*s^2-2*s^3;         C11=6*s-6*s^2;      C12=6-12*s;  C13=-12;
    C20=s-2*s^2+s^3;            C21=1-4*s+3*s^2;    C22=-4+6*s;  C23=6;
    C30=0;                      C31=0;              C32=0;       C33=0;
    C40=0;                      C41=0;              C42=0;       C43=0;

    G00 = H*(r-s)^3/6   +C10*r^3/6 +C20*r^2/2 +C30*r +C40;
    G01 = -H*(r-s)^2/2  +C11*r^3/6 +C21*r^2/2 +C31*r +C41;
    G02 = H*(r-s)       +C12*r^3/6 +C22*r^2/2 +C32*r +C42;
    G03 = -H + C13*r^3/6 +C23*r^2/2 +C33*r +C43;
    
    G20 = H*(r-s)       +C10*r      +C20;
    G21 = -H            +C11*r      +C21;  
    
    flag=h;
    if s==0 || s==1
        flag=flag/2;
    end
    G00=G00*flag;G01=G01*flag;G02=G02*flag;G20=G20*flag;
    G21=G21*flag; G03=G03*flag;
    end

    function [G00,G01,G02,G03,G20,G21]=MatrixOfGreen(j,k,h)

    % G(r,s)
    r=h*(j-1); 
    s=h*(k-1);
    H=(r>=s);
    
    %% Заделка - заделка
    C10=-1+3*s^2-2*s^3;         C11=6*s-6*s^2;      C12=6-12*s;  C13=-12;
    C20=s-2*s^2+s^3;            C21=1-4*s+3*s^2;    C22=-4+6*s;  C23=6;
    C30=0;                      C31=0;              C32=0;       C33=0;
    C40=0;                      C41=0;              C42=0;       C43=0;

    G00 = H*(r-s)^3/6    +C10*r^3/6 +C20*r^2/2 +C30*r +C40;
    G01 = -H*(r-s)^2/2   +C11*r^3/6 +C21*r^2/2 +C31*r +C41;
    G02 = H*(r-s)        +C12*r^3/6 +C22*r^2/2 +C32*r +C42;
    G03 = -H + C13*r^3/6 +C23*r^2/2 +C33*r +C43;
    
    G20 = H*(r-s)       +C10*r      +C20;
    G21 = -H            +C11*r      +C21;       
    end

    function [J]=MatrixJ(j,k,h)

    % J(r,s)
    r=h*(j-1); 
    s=h*(k-1);
 
    J = 6*r-12*s*r+6*s-4;
    flag=h;
    if s==0 || s==1
        flag=flag/2;
    end
    J=J*flag;
    end

end