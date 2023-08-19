function [G00Int,G01Int,G02Int,G00,G01,G02]=get_matrix(m)
    
    h = 1/m;
    
% Вычисление матриц функций
    for j=1:m+1     % проходим по строкам матриц
        for k=1:m+1 % проходим по столбцам матриц
            [G00Int(j,k),G01Int(j,k),G02Int(j,k)]=MatrixOfGreenIntegrate(j,k,h);
            [G00(j,k),G01(j,k),G02(j,k)]=MatrixOfGreen(j,k,h);
        end
    end
    
    function [G00,G01,G02]=MatrixOfGreenIntegrate(j,k,h)

    % G(r,s)
    r=h*(j-1); 
    s=h*(k-1);
    H=(r>=s);
    
    %% Заделка - заделка
    C10=-1+3*s^2-2*s^3;         C11=6*s-6*s^2;      C12=6-12*s;
    C20=s-2*s^2+s^3;            C21=1-4*s+3*s^2;    C22=-4+6*s;
    C30=0;                      C31=0;              C32=0;
    C40=0;                      C41=0;              C42=0;

    G00 = H*(r-s)^3/6   +C10*r^3/6 +C20*r^2/2 +C30*r +C40;
    G01 = -H*(r-s)^2/2   +C11*r^3/6 +C21*r^2/2 +C31*r +C41;
    G02 = H*(r-s)       +C12*r^3/6 +C22*r^2/2 +C32*r +C42; 
    
    flag=h;
    if s==0 || s==1
        flag=flag/2;
    end
    G00=G00*flag;G02=G02*flag;
    end

    function [G00,G01,G02]=MatrixOfGreen(j,k,h)

    % G(r,s)
    r=h*(j-1); 
    s=h*(k-1);
    H=(r>=s);
    
    %% Заделка - заделка
    C10=-1+3*s^2-2*s^3;         C11=6*s-6*s^2;      C12=6-12*s;
    C20=s-2*s^2+s^3;            C21=1-4*s+3*s^2;    C22=-4+6*s;
    C30=0;                      C31=0;              C32=0;
    C40=0;                      C41=0;              C42=0;

    G00 = H*(r-s)^3/6   +C10*r^3/6 +C20*r^2/2 +C30*r +C40;
    G01 = -H*(r-s)^2/2   +C11*r^3/6 +C21*r^2/2 +C31*r +C41;
    G02 = H*(r-s)       +C12*r^3/6 +C22*r^2/2 +C32*r +C42;
           
    end
end