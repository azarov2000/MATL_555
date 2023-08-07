function [G00Int,G01Int,G02Int,G10Int,G11Int,G12Int,...
          G00,G01,G02,G10,G11,G12]=get_matrix(m)
    
    h = 1/m;
    
% Вычисление матриц функций Грина
    for j=1:m+1     % проходим по строкам матриц Грина
        for k=1:m+1 % проходим по столбцам матриц Грина
            [G00Int(j,k),G01Int(j,k),G02Int(j,k),G10Int(j,k),G11Int(j,k),G12Int(j,k)]=MatrixOfGreenIntegrate(j,k,h);
            [G00(j,k),G01(j,k),G02(j,k),G10(j,k),G11(j,k),G12(j,k)]=MatrixOfGreen(j,k,h);
        end
    end
    

    function [G00,G01,G02,G10,G11,G12]=MatrixOfGreenIntegrate(j,k,h)

    % G(r,s)
    r=h*(j-1); 
    s=h*(k-1);
    H=(r>=s);
    
    %% Заделка - свободный край
    C10=-1;     C11=0;  C12=0;
    C20=s;      C21=1;  C22=0;
    C30=0;      C31=0;  C32=0;
    C40=0;      C41=0;  C42=0;

    G00 = H*(r-s)^3/6   +C10*r^3/6 +C20*r^2/2 +C30*r +C40;
    G01 = -H*(r-s)^2/2  +C11*r^3/6 +C21*r^2/2 +C31*r +C41;
    G02 = H*(r-s)       +C12*r^3/6 +C22*r^2/2 +C32*r +C42;
    
    G10 = H*(r-s)^2/2   +C10*r^2/2  +C20*r+ C30;
    G11 = -H*(r-s)      +C11*r^2/2  +C21*r+ C31;
    G12 = H             +C12*r^2/2  +C22*r+ C32;   
    
    flag=h;
    if s==0 || s==1
        flag=flag/2;
    end
    G00=G00*flag;G01=G01*flag;G02=G02*flag;G10=G10*flag;
    G11=G11*flag;G12=G12*flag;
    end

    function [G00,G01,G02,G10,G11,G12]=MatrixOfGreen(j,k,h)

    % G(r,s)
    r=h*(j-1); 
    s=h*(k-1);
    H=(r>=s);
    
    %% Заделка - свободный край
    C10=-1;     C11=0;  C12=0;
    C20=s;      C21=1;  C22=0;
    C30=0;      C31=0;  C32=0;
    C40=0;      C41=0;  C42=0;

    G00 = H*(r-s)^3/6   +C10*r^3/6 +C20*r^2/2 +C30*r +C40;
    G01 = -H*(r-s)^2/2  +C11*r^3/6 +C21*r^2/2 +C31*r +C41;
    G02 = H*(r-s)       +C12*r^3/6 +C22*r^2/2 +C32*r +C42;
    
    G10 = H*(r-s)^2/2   +C10*r^2/2  +C20*r+ C30;
    G11 = -H*(r-s)      +C11*r^2/2  +C21*r+ C31;
    G12 = H             +C12*r^2/2  +C22*r+ C32;    
    
    end
end