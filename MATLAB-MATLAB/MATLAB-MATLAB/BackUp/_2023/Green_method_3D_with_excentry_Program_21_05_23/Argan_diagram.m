function Argan_diagram (N,step_N,dat)

  
  flag = 0;    % флаг отсчёта 10 красных точек
  i = 1;       % счётчик цикла while
  mb_crit = 0; % начальное приближение для определения критической скорости
  while (flag<1)
      [A0,A1,A2] = maxtix_coeff_push_N(N,dat);
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
% ZEROFUN=@(x) interp1(N_vect(:),ReMax(:),x,'spline');
% Ncritical=fzero(ZEROFUN,mb_crit);

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
    function [A0,A1,A2] = maxtix_coeff_push_N(N,dat) 
    
        a_psi_0 = kron(dat.I,dat.E+2*dat.zeta_V*N*dat.S)-dat.Mom*kron(dat.G03Int,dat.S);
        a_X_0 = -12*dat.Mom*kron(dat.Fg,dat.E);

        b_psi_0 = 12*kron(dat.I(end,:),dat.E);
        b_X_0 = -12*dat.Mom*dat.Fg_ticks*dat.S+12*(dat.E+2*dat.zeta_V*N*dat.S)+6*dat.Mom*dat.S+dat.kappa_B*dat.E;

        A0 = [a_psi_0, a_X_0;
              b_psi_0, b_X_0];


        a_psi_1 = 2*dat.zeta_V*kron(dat.I,dat.E)-2*dat.zeta_e*kron(dat.G00Int,dat.E)-2*N*dat.beta*kron(dat.G02Int,dat.S);
        a_X_1 = -2*dat.zeta_e*kron(dat.F0,dat.E)-2*dat.beta*N*kron(dat.F2,dat.S);

        b_psi_1 = -2*dat.zeta_e*kron(dat.G30Int(end,:),dat.E)-2*N*dat.beta*kron(dat.T',dat.S);
        b_X_1 = -2*dat.zeta_e*dat.F0_ticks*dat.E-2*dat.beta*N*dat.F2_ticks*dat.S+24*dat.zeta_V*dat.E+2*dat.eta_B*dat.E;

        A1 = [a_psi_1, a_X_1;
              b_psi_1, b_X_1];

        a_psi_2 = kron(dat.G00Int,dat.E)-dat.beta*kron(dat.G02Int,dat.E);
        a_X_2 = kron(dat.F0,dat.E)-dat.beta*kron(dat.F2,dat.E);
        b_psi_2 = kron(dat.G30Int(end,:),dat.E)-dat.beta*kron(dat.T',dat.E);
        b_X_2 = dat.F0_ticks*dat.E-dat.beta*dat.F2_ticks*dat.E+dat.mu_B*dat.E; 

        A2 = [a_psi_2, a_X_2;
              b_psi_2, b_X_2];
    
    end 


end




