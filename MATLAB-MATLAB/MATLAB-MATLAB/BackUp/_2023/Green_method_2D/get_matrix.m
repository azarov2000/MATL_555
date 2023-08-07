function [Const_F1,Const_F2,Const_F1_ticks,Const_F2_ticks,Const_F3,...
          G00Int,G02Int,G30Int]=get_matrix(m)
    
    h = 1/m;
    z_dimless = linspace(0,1,1000);
    
% Вычисление векторов F1 и F2
    F1 = @(r,s)  ((r>s)*(r-s)^3/6   +(-1+3*s^2-2*s^3)*r^3/6 +(s-2*s^2+s^3)*r^2/2)*(6-12*s);
    F2 = @(r,s)  ((r>s)*(r-s)^3/6   +(-1+3*s^2-2*s^3)*r^3/6 +(s-2*s^2+s^3)*r^2/2)*(3*s^2-2*s^3);
    for i = 1:m+1
        for j = 1:length(z_dimless)
            f_1s{i}(j) = F1((i-1)*h,z_dimless(j));
            f_2s{i}(j) = F2((i-1)*h,z_dimless(j));
        end
        Const_F1(i) = trapz(z_dimless,f_1s{i}(:));
        Const_F2(i) = trapz(z_dimless,f_2s{i}(:));
    end 
    Const_F1=(Const_F1)';
    Const_F2=(Const_F2)';
    
% Вычисление строк F1''' и F2'''
    F1_ticks = @(r,s)  ((r>s)+(-1+3*s^2-2*s^3))*(6-12*s);
    F2_ticks = @(r,s)  ((r>s)+(-1+3*s^2-2*s^3))*(3*s^2-2*s^3);
    for j = 1:length(z_dimless)
        f_1s_ticks(j) = F1_ticks(1,z_dimless(j));
        f_2s_ticks(j) = F2_ticks(1,z_dimless(j));
    end
    Const_F1_ticks = trapz(z_dimless,f_1s_ticks(:));
    Const_F2_ticks = trapz(z_dimless,f_2s_ticks(:));
    
% Вычисление вспомогительной функции для матрицы G32
    F3 = @(s) (6-12*s);
    for i = 1:m+1
        Const_F3(i) = F3((i-1)*h);  
    end
% Вычисление матриц функций Грина
    for j=1:m+1     % проходим по строкам матриц Грина
        for k=1:m+1 % проходим по столбцам матриц Грина
            [G00Int(j,k),G02Int(j,k),G30Int(j,k)]=MatrixOfGreenIntegrate(j,k,h);
        end
    end
    
    function [G00,G02,G30]=MatrixOfGreenIntegrate(j,k,h)

    % G(r,s)

    r=h*(j-1); 
    s=h*(k-1);
    H=(r>s);
    %% Заделка - заделка
    C10=-1+3*s^2-2*s^3;     C12=6-12*s;
    C20=s-2*s^2+s^3;        C22=-4+6*s;
    C30=0;                  C32=0;
    C40=0;                  C42=0;

    G00 = H*(r-s)^3/6   +C10*r^3/6 +C20*r^2/2 +C30*r +C40;
    G02 = H*(r-s)       +C12*r^3/6 +C22*r^2/2 +C32*r +C42;
    G30 = H + C10;
    
    flag=h;
    if s==0 || s==1
        flag=flag/2;
    end
    G00=G00*flag;G02=G02*flag;G30=G30*flag;
    end
end