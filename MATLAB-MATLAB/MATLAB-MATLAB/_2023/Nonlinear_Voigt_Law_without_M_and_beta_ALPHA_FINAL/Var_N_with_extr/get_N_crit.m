function Ncritical = get_N_crit (N,step_N,d)
  
  flag = 0;    % флаг отсчёта 10 красных точек
  i = 1;       % счётчик цикла while
  mb_crit = 0; % начальное приближение для определения критической скорости
  while (flag<3)
      [A0,A1,A2] = maxtix_for_Argan(N,d);
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
Ncritical=fzero(ZEROFUN,mb_crit);

% Построение диаграммы Аргана
% figure; 
% subplot(211)
% ff = gca; 
% ff.FontName = 'Times New Roman';
% ff.FontSize = 20;
% hold on; box on; grid on; 
% for j = 1:length(N_vect)
%     marker_type = '.k';
%     marker_size = 17;
%     if max(real(EigenValues{j}))>0
%        marker_type = '.r';
%        marker_size = 20;
%     end
%     plot (N_vect(j), real(EigenValues{j}),marker_type,'MarkerSize',marker_size)
% end
% xlabel('N')
% ylabel('Re(\lambda)')
% 
% subplot(212)
% ff = gca; 
% ff.FontName = 'Times New Roman';
% ff.FontSize = 20;
% hold on; box on; grid on; 
% for j = 1:length(N_vect)
%     marker_type = '.k';
%     marker_size = 17;
%     if max(real(EigenValues{j}))>0
%        marker_type = '.r';
%        marker_size = 20;
%     end
%     plot (N_vect(j), imag(EigenValues{j}),marker_type,'MarkerSize',marker_size)
% end
% xlabel('N')
% ylabel('Im(\lambda)')
% 
% figure;
% ff = gca; 
% ff.FontName = 'Times New Roman';
% ff.FontSize = 20;
% hold on; box on; grid on; 
% for j = 1:length(N_vect)
%     marker_type = '.k';
%     marker_size = 17;
%     if max(real(EigenValues{j}))>0
%        marker_type = '.r';
%        marker_size = 20;
%     end
%     plot (real(EigenValues{j}), imag(EigenValues{j}),marker_type,'MarkerSize',marker_size)
% end
% xlabel('Re(\lambda)')
% ylabel('Im(\lambda)')

 disp(['Критическая скорость равна ',num2str(Ncritical)]);


% Функция возвращает матрицы коэффициентов при заданном N
    function [A0,A1,A2] = maxtix_for_Argan(N,d)

        A0 = kron(d.I,(d.E-2*d.zeta_V*N*d.R));

        A1 = 2*d.zeta_V*kron(d.I,d.E)+2*d.zeta_e*kron(d.G00Int,d.E);

        A2 = d.mu_R*kron(d.G00Int,d.E);
    
    end 


end




