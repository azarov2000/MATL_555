function Argan_diagram (N,step_N)

    global I E S zeta_V Mom G03Int
    global Const_Fg IE Const_Fg_ticks kappa_B 
    global zeta_e G00Int beta G02Int Const_F0 Const_F2
    global G30E Const_T Const_F0_ticks Const_F2_ticks
    global eta_B mu_B
  
  flag = 0;    % флаг отсчёта 10 красных точек
  i = 1;       % счётчик цикла while
  mb_crit = 0; % начальное приближение для определения критической скорости
  while (flag<1)
      [A0,A1,A2] = maxtix_for_Argan(N);
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
figure; 
subplot(211)
ff = gca; 
ff.FontName = 'Times New Roman';
ff.FontSize = 20;
hold on; box on; grid on; 
for j = 1:length(N_vect)
    marker_type = '.k';
    marker_size = 17;
    if max(real(EigenValues{j}))>0
       marker_type = '.r';
       marker_size = 20;
    end
    plot (N_vect(j), real(EigenValues{j}),marker_type,'MarkerSize',marker_size)
end
xlabel('N')
ylabel('Re(\lambda)')

subplot(212)
ff = gca; 
ff.FontName = 'Times New Roman';
ff.FontSize = 20;
hold on; box on; grid on; 
for j = 1:length(N_vect)
    marker_type = '.k';
    marker_size = 17;
    if max(real(EigenValues{j}))>0
       marker_type = '.r';
       marker_size = 20;
    end
    plot (N_vect(j), imag(EigenValues{j}),marker_type,'MarkerSize',marker_size)
end
xlabel('N')
ylabel('Im(\lambda)')

figure;
ff = gca; 
ff.FontName = 'Times New Roman';
ff.FontSize = 20;
hold on; box on; grid on; 
for j = 1:length(N_vect)
    marker_type = '.k';
    marker_size = 17;
    if max(real(EigenValues{j}))>0
       marker_type = '.r';
       marker_size = 20;
    end
    plot (real(EigenValues{j}), imag(EigenValues{j}),marker_type,'MarkerSize',marker_size)
end
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')

 disp(['Критическая скорость равна ',num2str(Ncritical)]);


% Функция возвращает матрицы коэффициентов при заданном N
    function [A0,A1,A2] = maxtix_for_Argan(N) 
    
A0 = [kron(I,E+2*zeta_V*N*S)-Mom*kron(G03Int,S),   -12*Mom*kron(Const_Fg,E);
      12*Mom*IE(end-1:end,:),                            -12*Mom*Const_Fg_ticks*S+12*(E+2*zeta_V*N*S)+6*Mom*S+kappa_B*E];

A1 = [2*zeta_V*kron(I,E)-2*zeta_e*kron(G00Int,E)-2*N*beta*kron(G02Int,S),   -2*zeta_e*kron(Const_F0,E)-2*beta*N*kron(Const_F2,S);
      -2*zeta_e*G30E(end-1:end,:)-2*N*beta*kron(Const_T',S),                -2*zeta_e*Const_F0_ticks*E-2*beta*N*Const_F2_ticks*S+24*zeta_V*E+2*eta_B*E];

A2 = [kron(G00Int,E)-beta*kron(G02Int,E),       kron(Const_F0,E)-beta*kron(Const_F2,E);
      G30E(end-1:end,:)-beta*kron(Const_T',E),  Const_F0_ticks*E-beta*Const_F2_ticks*E+mu_B*E];
    
    end 


end




