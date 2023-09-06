function Ncritical = get_N_crit (dN,step_N,d)
  
  flag = 0;    % флаг отсчёта 10 красных точек
  i = 1;       % счётчик цикла while
  mb_crit = 0; % начальное приближение для определения критической скорости
  while (flag<2)
      d.N = dN;
      [A0,A1,A2] = Filling_matrix_A(d);
      EigenValues{i} = polyeig(A0,A1,A2);
      if max(real(EigenValues{i})>0)
          flag = flag+1;
      end
      if max(real(EigenValues{i})>0) && mb_crit==0
          mb_crit = dN;
      end 
      N_vect(i) = dN;
      dN = dN+step_N;
      i=i+1;
  end
% Определение критической скорости вращения
for j=1:length(N_vect)
    ReMax(j)=max(real(EigenValues{j}));
end
% ZEROFUN=@(x) interp1(N_vect(:),ReMax(:),x,'spline');
ZEROFUN=@(x) interp1(N_vect(:),ReMax(:),x);
Ncritical=fzero(ZEROFUN,mb_crit);

disp(['Критическая скорость равна ',num2str(Ncritical)]);
 


end




