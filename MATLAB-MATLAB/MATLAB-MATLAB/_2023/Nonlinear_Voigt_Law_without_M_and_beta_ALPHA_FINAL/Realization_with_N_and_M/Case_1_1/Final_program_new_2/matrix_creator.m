function [mc]=matrix_creator(m,zeta_e,zeta_V,N_z,M_z)
    
h = 1/m;
    
G00Int = zeros(m+1);
G01Int = zeros(m+1);
G02Int = zeros(m+1);
G00 = zeros(m+1);
G01 = zeros(m+1);
G02 = zeros(m+1);

% Вычисление матриц функций Грина
for f=1:m+1     % проходим по строкам матриц
    for k=1:m+1 % проходим по столбцам матриц
        [G00Int(f,k),G01Int(f,k),G02Int(f,k)]=MatrixOfGreenIntegrate(f,k,h);
        [G00(f,k),G01(f,k),G02(f,k)]=MatrixOfGreen(f,k,h);
    end
end

% Создание дополнительных матриц
I=eye(m+1);               
Z=zeros(m+1);
mc = edging(G00Int,G01Int,G02Int,G00,G01,G02,I,Z,m);
mc.E = eye(2);  
mc.R = [0 -1; 1 0];

mc.m = m;
mc.zeta_e = zeta_e;
mc.zeta_V = zeta_V;
mc.N_z = N_z;
mc.M_z = M_z;
mc.M = get_matr_coeff_M(m);


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

function [mc] = edging(G00Int,G01Int,G02Int,G00,G01,G02,I,Z,m)

    Lim = m;
    mc.G00Int = G00Int(2:Lim,2:Lim);
    mc.G01Int = G01Int(2:Lim,2:Lim);
    mc.G02Int = G02Int(2:Lim,2:Lim);
    
    mc.G00 = G00(2:Lim,2:Lim);
    mc.G01 = G01(2:Lim,2:Lim);
    mc.G02 = G02(2:Lim,2:Lim);
    
    mc.I = I(2:Lim,2:Lim);
    mc.Z = Z(2:Lim,2:Lim);

end

function M = get_matr_coeff_M(m)
    h = 1/m;
    M = zeros(m-1);
    
    for i=1:length(M)
        for j=1:length(M)
            
            if i==j
                M(i,j) = -2;
            end
            if i == j-1 || j == i-1 
                M(i,j) = 1;
            end 
            
        end 
    end
    M = (1/h)^2*M;
    M = kron(M,eye(2));
end


end