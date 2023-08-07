function Argan_diagram (N,step_N)

    global I E S zeta_V Mom 
    global  kappa_B 
    global zeta_e G00Int beta G02Int
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
delta = 1e-8;

a_ksi_0 = kron(I,(E+2*zeta_V*N*S))+kappa_B*kron(G00*EC,E);
a_theta_0 = Mom*(h^-1)*kron(G00*EC,E)-Mom*kron(G02Int,E);

b_ksi_0 = -kappa_B*kron(G10*EC,S);
b_theta_0 = kron(I,(E+2*zeta_V*N*S))-Mom*(h^-1)*kron(G10*EC,S)+Mom*kron(G12Int,S);

% new
b_ksi_0 = b_ksi_0 + delta*ones(length(b_ksi_0));
b_theta_0 = b_theta_0 + delta*ones(length(b_theta_0));

A0 = [a_ksi_0, a_theta_0;
      b_ksi_0, b_theta_0];

a_ksi_1 = 2*zeta_V*kron(I,E)+2*zeta_e*kron(G00Int,E)+2*eta_B*kron(G00*EC,E);
a_theta_1 = -2*beta*N*kron(G01Int,E);

b_ksi_1 = -2*zeta_e*kron(G10Int,S)-2*eta_B*kron(G10*EC,S);
b_theta_1 = 2*zeta_V*kron(I,E)+2*beta*N*kron(G11Int,S);

% new
b_ksi_1 = b_ksi_1+delta*ones(length(b_ksi_1));
b_theta_1 = b_theta_1+delta*ones(length(b_theta_1));

A1 = [a_ksi_1, a_theta_1;
      b_ksi_1, b_theta_1];

  
a_ksi_2 = kron(G00Int,E)+mu_B*kron(G00*EC,E);
a_theta_2 = beta*kron(G01Int,S);

b_ksi_2 = -kron(G10Int,S)-mu_B*kron(G10*EC,S);
b_theta_2 = beta*kron(G11Int,E);

% new
b_ksi_2 = b_ksi_2+delta*ones(length(b_ksi_2));
b_theta_2 = b_theta_2+delta*ones(length(b_theta_2));

A2 = [a_ksi_2, a_theta_2;
      b_ksi_2, b_theta_2];
    
    end 


end




