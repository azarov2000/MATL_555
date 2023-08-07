function [Const_F00,Const_F02,Const_Fg0,Const_F10,Const_F12,Const_Fg1,...
          Const_F10_ticks,Const_F12_ticks,Const_Fg1_ticks,...
          Const_T,Const_P,...
          G00Int,G10Int,G01Int,G11Int,G02Int,G12Int,G30Int]=get_matrix_new(m)
    


      
   
      
    h = 1/m;
    z_dimless = linspace(0,1,10000);

% Вычисление векторов F00, F02, Fg0, F10, F12, Fg1
    F00 = @(r,s)  ((r>s)*(r-s)^3/6   +(-1+3*s^2-2*s^3)*r^3/6 +(s-2*s^2+s^3)*r^2/2)*(3*s^2-2*s^3);
    F02 = @(r,s)  ((r>s)*(r-s)^3/6   +(-1+3*s^2-2*s^3)*r^3/6 +(s-2*s^2+s^3)*r^2/2)*(6-12*s);
    Fg0 = @(r,s)  ((r>s)*(r-s)^3/6   +(-1+3*s^2-2*s^3)*r^3/6 +(s-2*s^2+s^3)*r^2/2);

    F10 = @(r,s)  ((r>s)*(r-s)^2/2   +(-1+3*s^2-2*s^3)*r^2/2 +(s-2*s^2+s^3)*r)*(3*s^2-2*s^3);
    F12 = @(r,s)  ((r>s)*(r-s)^2/2   +(-1+3*s^2-2*s^3)*r^2/2 +(s-2*s^2+s^3)*r)*(6-12*s);
    Fg1 = @(r,s)  ((r>s)*(r-s)^2/2   +(-1+3*s^2-2*s^3)*r^2/2 +(s-2*s^2+s^3)*r);
    
    
    for i = 1:m+1
        for j = 1:length(z_dimless)
            f_00s{i}(j) = F00((i-1)*h,z_dimless(j));
            f_02s{i}(j) = F02((i-1)*h,z_dimless(j));
            f_g0s{i}(j) = Fg0((i-1)*h,z_dimless(j));
            
            f_10s{i}(j) = F10((i-1)*h,z_dimless(j));
            f_12s{i}(j) = F12((i-1)*h,z_dimless(j));
            f_g1s{i}(j) = Fg1((i-1)*h,z_dimless(j));      
        end
        Const_F00(i) = trapz(z_dimless,f_00s{i}(:));
        Const_F02(i) = trapz(z_dimless,f_02s{i}(:));
        Const_Fg0(i) = trapz(z_dimless,f_g0s{i}(:));

        Const_F10(i) = trapz(z_dimless,f_10s{i}(:));
        Const_F12(i) = trapz(z_dimless,f_12s{i}(:));
        Const_Fg1(i) = trapz(z_dimless,f_g1s{i}(:));
    end 
    Const_F00=(Const_F00)';
    Const_F02=(Const_F02)';
    Const_Fg0=(Const_Fg0)';

    Const_F10=(Const_F10)';
    Const_F12=(Const_F12)';
    Const_Fg1=(Const_Fg1)';
    
    
% Вычисление строк F10'', F12'', Fg1''

    F10_ticks = @(r,s)  ((r>s)+(-1+3*s^2-2*s^3))*(3*s^2-2*s^3);
    F12_ticks = @(r,s)  ((r>s)+(-1+3*s^2-2*s^3))*(6-12*s);
    Fg1_ticks = @(r,s)  ((r>s)+(-1+3*s^2-2*s^3));
    
    for j = 1:length(z_dimless)
        f_10s_ticks(j) = F10_ticks(1,z_dimless(j));
        f_12s_ticks(j) = F12_ticks(1,z_dimless(j));
        f_g1s_ticks(j) = Fg1_ticks(1,z_dimless(j));
    end
    Const_F10_ticks = trapz(z_dimless,f_10s_ticks(:));
    Const_F12_ticks = trapz(z_dimless,f_12s_ticks(:));
    Const_Fg1_ticks = trapz(z_dimless,f_g1s_ticks(:));
    
% Вычисление векторов T и P
    T = @(s) (6-12*s);
    P = @(s) (6*s-6*s^2);
    for i = 1:m+1
        Const_T(i) = T((i-1)*h);
        Const_P(i) = P((i-1)*h);
    end
    Const_T = Const_T';
    Const_P = Const_P';
    

% Вычисление матриц функций Грина
    for j=1:m+1     % проходим по строкам матриц Грина
        for k=1:m+1 % проходим по столбцам матриц Грина
            [G00Int(j,k),G10Int(j,k),G01Int(j,k),G11Int(j,k),G02Int(j,k),G12Int(j,k),G30Int(j,k)]=...
                MatrixOfGreenIntegrate(j,k,h);
        end
    end
    
    function [G00,G10,G01,G11,G02,G12,G30]=MatrixOfGreenIntegrate(j,k,h)

    % G(r,s)

    r=h*(j-1); 
    s=h*(k-1);
    H=(r>=s);
    %% Заделка - заделка
    C10=-1+3*s^2-2*s^3;  C11=6*s-6*s^2;           C12=6-12*s; C13=-12;
    C20=s-2*s^2+s^3;     C21=1-4*s+3*s^2;         C22=-4+6*s; C23=6; 
    C30=0;               C31=0;                   C32=0;      C33=0;
    C40=0;               C41=0;                   C42=0;      C43=0;

    G00 = H*(r-s)^3/6   +C10*r^3/6 +C20*r^2/2 +C30*r +C40;
    G10 = H*(r-s)^2/2   +C10*r^2/2 +C20*r     +C30;
    G01 = -H*(r-s)^2/2  +C11*r^3/6 +C21*r^2/2 +C31*r +C41;
    G11 = -H*(r-s)      +C11*r^2/2 +C21*r     +C31;
    
    G02 = H*(r-s)       +C12*r^3/6 +C22*r^2/2 +C32*r +C42;
    G12 = H +C12*r^2/2 +C22*r +C32;
    
    G30 = H + C10;
    
    flag=h;
    if s==0 || s==1
        flag=flag/2;
    end
    G00=G00*flag;G10=G10*flag;G01=G01*flag;G11=G11*flag;
    G02=G02*flag;G12=G12*flag;G30=G30*flag;
    end
end