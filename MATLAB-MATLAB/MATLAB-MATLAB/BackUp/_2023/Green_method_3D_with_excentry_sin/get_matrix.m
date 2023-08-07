function [F_00, F_02, F_03, F_4, F_30C, F_32C, F_33C, T,...
          G00Int,G02Int,G03Int,G30Int]=get_matrix(m)
    
h = 1/m;
z_dimless = linspace(0,1,10000);
    
% Вычисление векторов
F_00_f = @(r,s)  ((r>=s)*(r-s)^3/6   +(-1+3*s^2-2*s^3)*r^3/6 +(s-2*s^2+s^3)*r^2/2)*0.5*(1-cos(pi*s));
F_02_f = @(r,s)  ((r>=s)*(r-s)^3/6   +(-1+3*s^2-2*s^3)*r^3/6 +(s-2*s^2+s^3)*r^2/2)*(pi^2/2)*cos(pi*s);
F_03_f = @(r,s)  ((r>=s)*(r-s)^3/6   +(-1+3*s^2-2*s^3)*r^3/6 +(s-2*s^2+s^3)*r^2/2)*(-pi^3/2)*sin(pi*s);
F4_f = @(r) -(pi^4/2)*cos(pi*r);
T_f = @(r) (6-12*r);
for i = 1:m+1
    for j = 1:length(z_dimless)
        f_00{i}(j) = F_00_f((i-1)*h,z_dimless(j));
        f_02{i}(j) = F_02_f((i-1)*h,z_dimless(j));
        f_03{i}(j) = F_03_f((i-1)*h,z_dimless(j));
    end
    F_00(i) = trapz(z_dimless,f_00{i}(:));
    F_02(i) = trapz(z_dimless,f_02{i}(:));
    F_03(i) = trapz(z_dimless,f_03{i}(:));
    F_4(i) = F4_f((i-1)*h);
    T(i) = T_f((i-1)*h);
end
F_00 = F_00'; F_02 = F_02'; F_03 = F_03'; F_4 = F_4'; T = T';

% Вычисление констант
F_30C_f = @(r,s) ((r>=s)*(-1+3*s^2-2*s^3))*0.5*(1-cos(pi*s));
F_32C_f = @(r,s) ((r>=s)*(-1+3*s^2-2*s^3))*(pi^2/2)*cos(pi*s);
F_33C_f = @(r,s) ((r>=s)*(-1+3*s^2-2*s^3))*(-pi^3/2)*sin(pi*s);


for j = 1:length(z_dimless)
    f_30C(j) = F_30C_f(1,z_dimless(j));
    f_32C(j) = F_32C_f(1,z_dimless(j));
    f_33C(j) = F_33C_f(1,z_dimless(j));
end
F_30C = trapz(z_dimless,f_30C(:));
F_32C = trapz(z_dimless,f_32C(:));
F_33C = trapz(z_dimless,f_33C(:));

F_30C = F_30C';F_32C = F_32C'; F_33C = F_33C';


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
    H=(r>=s);
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