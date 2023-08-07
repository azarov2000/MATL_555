close all
clc
clear
%%
% Parameters
d = 20*10^-3;     %[m]
l = 0.7;            %[m]
ro = 7800;          %[kg/m^3]

N = 23;              %[-]
eps_d = d/l;        %[-]
beta = (eps_d^2)/16;%[-]
zeta_e = 0.025;     %[-]
zeta_V = 0.005;     %[-]


rootS = solver_for_alpha(1,1,20);
alpha1 = rootS(1);
E = eye(2); S=[0 1;-1 0];
A0 = alpha1^4 * (E+2*zeta_V*N*S);
A1 = 2*(alpha1^4*zeta_V+zeta_e)*E;
A2 = E;



%% Поиск скритической скорости

    N = 20;          % начальная скорость
    step_N = 2;  % шаг по скоростям
    flag = 0;       % флаг отсчёта 3 красных точек
    i = 1;          % счётчик цикла while
    mb_crit = 0;    % начальное приближение для определения критической скорости
    while (flag<3)
      [A0,A1,A2]=test_matrix(N);
      EigenValues{i} = polyeig(A0,A1,A2);
      if max(real(EigenValues{i})>0)
          flag = flag+1;
      end
      if max(real(EigenValues{i})>0) && mb_crit==0
          mb_crit = N;
      end 
      N_vect(i) = N;
      N = N+step_N;
      i=i+1;
    end
    
    % Определение критической скорости вращения
    for j=1:length(N_vect)
        ReMax(j)=max(real(EigenValues{j}));
    end
    ZEROFUN=@(x) interp1(N_vect(:),ReMax(:),x,'spline');
    N_critical=fzero(ZEROFUN,mb_crit)




