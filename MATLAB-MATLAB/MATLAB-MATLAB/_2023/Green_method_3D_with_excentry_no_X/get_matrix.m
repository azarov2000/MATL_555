function [Const_F0,Const_F2,Const_Fg,Const_Fg30,Const_F0_ticks,Const_F2_ticks,Const_Fg_ticks,Const_T,...
          G00Int,G02Int,G03Int,G30Int]=get_matrix(m)
    
    h = 1/m;
    z_dimless = linspace(0,1,100000);
    
% Вычисление векторов F1 и F2
    F0 = @(r,s)  ((r>s)*(r-s)^3/6   +(-1+3*s^2-2*s^3)*r^3/6 +(s-2*s^2+s^3)*r^2/2)*(3*s^2-2*s^3);
    F2 = @(r,s)  ((r>s)*(r-s)^3/6   +(-1+3*s^2-2*s^3)*r^3/6 +(s-2*s^2+s^3)*r^2/2)*(6-12*s);
    Fg = @(r,s)  ((r>s)*(r-s)^3/6   +(-1+3*s^2-2*s^3)*r^3/6 +(s-2*s^2+s^3)*r^2/2);
    for i = 1:m+1
        for j = 1:length(z_dimless)
            f_0s{i}(j) = F0((i-1)*h,z_dimless(j));
            f_2s{i}(j) = F2((i-1)*h,z_dimless(j));
            f_gs{i}(j) = Fg((i-1)*h,z_dimless(j));
        end
        Const_F0(i) = trapz(z_dimless,f_0s{i}(:));
        Const_F2(i) = trapz(z_dimless,f_2s{i}(:));
        Const_Fg(i) = trapz(z_dimless,f_gs{i}(:));
    end 
    Const_F0=(Const_F0)';
    Const_F2=(Const_F2)';
    Const_Fg=(Const_Fg)';
    
% Вычисление строк F1''' и F2'''
    F0_ticks = @(r,s)  ((r>s)+(-1+3*s^2-2*s^3))*(3*s^2-2*s^3);
    F2_ticks = @(r,s)  ((r>s)+(-1+3*s^2-2*s^3))*(6-12*s);
    Fg_ticks = @(r,s)  ((r>s)+(-1+3*s^2-2*s^3));
    Fg30 = @(r,s) (r>s)*(-1+3*s^2-2*s^3);
    for j = 1:length(z_dimless)
        f_0s_ticks(j) = F0_ticks(1,z_dimless(j));
        f_2s_ticks(j) = F2_ticks(1,z_dimless(j));
        f_gs_ticks(j) = Fg_ticks(1,z_dimless(j));
        f_g30s(j) = Fg30(1,z_dimless(j));
    end
    Const_F0_ticks = trapz(z_dimless,f_0s_ticks(:));
    Const_F2_ticks = trapz(z_dimless,f_2s_ticks(:));
    Const_Fg_ticks = trapz(z_dimless,f_gs_ticks(:));
    Const_Fg30 = trapz(z_dimless,f_g30s(:));
    
% Вычисление вспомогительной функции для матрицы G32
    T = @(s) (6-12*s);
    for i = 1:m+1
        Const_T(i) = T((i-1)*h);  
    end
    Const_T = Const_T';
% Вычисление матриц функций Грина
    for j=1:m+1     % проходим по строкам матриц Грина
        for k=1:m+1 % проходим по столбцам матриц Грина
            [G00Int(j,k),G02Int(j,k),G03Int(j,k),G30Int(j,k)]=MatrixOfGreenIntegrate(j,k,h);
        end
    end
    
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
end