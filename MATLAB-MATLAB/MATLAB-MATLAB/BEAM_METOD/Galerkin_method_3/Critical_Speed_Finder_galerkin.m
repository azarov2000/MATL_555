function N_critical = Critical_Speed_Finder_galerkin(type_BC,N_start,N_step,n,beta,beta_R,mu_R,zeta_V,zeta_e,zeta_Dksi,zeta_Dte,zeta_C)
    N = N_start;          % начальная скорость
    step_N = N_step;      % шаг по скоростям
    flag = 0;             % флаг отсчёта 3 красных точек
    i = 1;                % счётчик цикла while
    mb_crit = 0;          % начальное приближение для определения критической скорости
    roots_alpha = solver_for_alpha(type_BC,1,50);
    while (flag<3)
      [A0,A1,A2] = get_matrix(type_BC,N,n,beta,beta_R,mu_R,zeta_V,zeta_e,zeta_Dksi,zeta_Dte,zeta_C,roots_alpha);
      EigenValues{i} = polyeig(A0,A1,A2);
      max(real(EigenValues{i})
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
    N_critical=fzero(ZEROFUN,mb_crit);
end