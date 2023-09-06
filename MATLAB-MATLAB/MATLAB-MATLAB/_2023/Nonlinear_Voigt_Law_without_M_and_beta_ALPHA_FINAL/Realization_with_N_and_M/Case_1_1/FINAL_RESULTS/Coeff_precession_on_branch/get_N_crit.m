function [Ncritical,freq_biff] = get_N_crit (dN,step_N,d)
  
  flag = 0;    % флаг отсчёта 10 красных точек
  i = 1;       % счётчик цикла while
  mb_crit = 0; % начальное приближение для определения критической скорости
  while (flag<3)
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
ZEROFUN=@(x) interp1(N_vect(:),ReMax(:),x,'spline');
Ncritical=fzero(ZEROFUN,mb_crit);

disp(['Критическая скорость равна ',num2str(Ncritical)]);

% Частота при бифуркации
tmp_mc = d;
tmp_mc.N = Ncritical;
[A0,A1,A2] = Filling_matrix_A(tmp_mc);
EigenVal = polyeig(A0,A1,A2);


% index = find(max(real(EigenVal))>-10e6 & max(real(EigenVal))<10e6);

index = find(max(real(EigenVal))==real(EigenVal));
freq_biff = imag(EigenVal(index(1)));
end




