clc
clear
close all

m = 6;          % кол-во участков разбиения
h = 1/m;        % шаг разбиения
N = 0;          % [-] Скорость вращения
eta_V = 0.005;  % [-] Коэфф. внешнего трения
 
% Матрицы 
I=eye(m+1);
E = eye(2);
R = [0 -1; 1 0];
G03Int = zeros(m+1);
for j=1:m+1     % проходим по строкам матрицы
    for k=1:m+1 % проходим по столбцам матрицы
        [G03Int(j,k)]=MatrixOfGreenIntegrate(j,k,h);
    end
end

%% Вычисление определителей при разном крутящем моменте

vector_M = linspace(-10,10,1000);
for i=1:length(vector_M)
    [DET(i)] = get_det(I,E,eta_V,N,R,G03Int,vector_M(i));
end



%% График
figure;
grid on; box on; hold on;
for i=1:length(vector_M)
    if DET(i)<=0
       color_marker = '.r';
    else
        color_marker = '.k';
    end
    plot(vector_M(i),DET(i),color_marker,'MarkerSize', 14)
end
ff = gca;
ff.FontSize = 20;
xlabel('M')
ylabel('det [ A_{ 0} ]')



function [G03]=MatrixOfGreenIntegrate(j,k,h)
    % G(r,s)
    r=h*(j-1); 
    s=h*(k-1);
    H=(r>=s);
    %% Заделка - заделка
    C10=-1+3*s^2-2*s^3;         C11=6*s-6*s^2;      C12=6-12*s;  C13=-12;
    C20=s-2*s^2+s^3;            C21=1-4*s+3*s^2;    C22=-4+6*s;  C23=6;
    C30=0;                      C31=0;              C32=0;       C33=0;
    C40=0;                      C41=0;              C42=0;       C43=0;
    
    G03 = -H + C13*r^3/6 +C23*r^2/2 +C33*r +C43;

    flag=h;
    if s==0 || s==1
        flag=flag/2;
    end
    G03=G03*flag;
end

function [DET]=get_det(I,E,eta_V,N,R,G03Int,M)
    A = kron(I,E)-2*eta_V*N*kron(I,R)+M*kron(G03Int,R);
    DET = det(A);
end
