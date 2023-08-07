function Argan_diagram (N,step_N)

    global I E zeta_V S G02Int Mom Fg0 Z G12Int Fg1
    global Fg1_ticks kappa_B G00Int zeta_e beta F00 F10 F02 F12
    global eta_B F12_ticks G10Int G30Int G11Int G01Int F10_ticks
    global mu_B P T
    
  flag = 0;    % флаг отсчёта 10 красных точек
  i = 1;       % счётчик цикла while
  mb_crit = 0; % начальное приближение для определения критической скорости
  while (flag<10)
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
        
A0 = [kron(I,E+2*zeta_V*N*S), -Mom*kron(G02Int,E), -12*Mom*kron(Fg0,S);
      kron(Z,E), kron(I,(E+2*zeta_V*N*S))+Mom*kron(G12Int,S), -12*Mom*kron(Fg1,E);
      kron(Z(end,:),E), -Mom*kron(T',E), -12*Mom*Fg1_ticks*S+12*(E+2*zeta_V*N*S)+6*Mom*S+kappa_B*E];



A1 = [2*zeta_V*kron(I,E)+2*zeta_e*kron(G00Int,E), -2*beta*N*kron(G01Int,E), 2*zeta_e*kron(F00,E)-2*beta*N*kron(F02,S);
      -2*zeta_e*kron(G10Int,S), 2*zeta_V*kron(I,E)+2*beta*N*kron(G11Int,S), -2*zeta_e*kron(F10,S)-2*beta*N*kron(F12,E);
      2*zeta_e*kron(G30Int(end,:),E), -2*beta*N*kron(P',E), 2*zeta_e*F10_ticks*E-2*beta*N*F12_ticks*S+24*zeta_V*E+2*eta_B*E];



A2 = [kron(G00Int,E), beta*kron(G01Int,S), kron(F00,E)-beta*kron(F02,E);
      -kron(G10Int,S), beta*kron(G11Int,E), -kron(F10,S)+beta*kron(F12,S);
      kron(G30Int(end,:),E), beta*kron(P',S), mu_B*E+F10_ticks*E-beta*F12_ticks*E];
    
    end 


end




