function Ncritical=Argan_diagram (data,N,step_N)

%     global I E S zeta_V Mom
%     global kappa_B 
%     global zeta_e G00Int beta G02Int
%     global eta_B mu_B
%     global EC G00
%     global G12Int G10 G01Int G11Int G10Int G01 G11
  flag = 0;    % флаг отсчёта 10 красных точек
  i = 1;       % счётчик цикла while
  mb_crit = 0; % начальное приближение для определения критической скорости
  while (flag<20)
      [A0,A1,A2] = maxtix_for_Argan(data,N);
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

disp(['Критическая скорость равна ',num2str(Ncritical)]);

% Функция возвращает матрицы коэффициентов при заданном N
function [A0,A1,A2] = maxtix_for_Argan(data,N) 
    
a_ksi_0 = kron(data.I,(data.E+2*data.zeta_V*N*data.S))+data.kappa_B*kron(data.G00*data.EC,data.E);
a_theta_0 = data.Mom*kron(data.G01*data.EC,data.E)-data.Mom*kron(data.G02Int,data.E);

b_ksi_0 = -data.kappa_B*kron(data.G10*data.EC,data.S);
b_theta_0 = kron(data.I,(data.E+2*data.zeta_V*N*data.S))-data.Mom*kron(data.G11*data.EC,data.S)+data.Mom*kron(data.G12Int,data.S);

A0 = [a_ksi_0, a_theta_0;
      b_ksi_0, b_theta_0];

a_ksi_1 = 2*data.zeta_V*kron(data.I,data.E)+2*data.zeta_e*kron(data.G00Int,data.E)+2*data.eta_B*kron(data.G00*data.EC,data.E);
a_theta_1 = 2*data.beta*N*kron(data.G00*data.EC,data.E)-2*data.beta*N*kron(data.G01Int,data.E);

b_ksi_1 = -2*data.zeta_e*kron(data.G10Int,data.S)-2*data.eta_B*kron(data.G10*data.EC,data.S);
b_theta_1 = 2*data.zeta_V*kron(data.I,data.E)-2*data.beta*N*kron(data.G10*data.EC,data.S)+2*data.beta*N*kron(data.G11Int,data.S);

A1 = [a_ksi_1, a_theta_1;
      b_ksi_1, b_theta_1];

a_ksi_2 = kron(data.G00Int,data.E)+data.mu_B*kron(data.G00*data.EC,data.E);
a_theta_2 = -data.beta*kron(data.G00*data.EC,data.S)+data.beta*kron(data.G01Int,data.S);

b_ksi_2 = -kron(data.G10Int,data.S)-data.mu_B*kron(data.G10*data.EC,data.S);
b_theta_2 = -data.beta*kron(data.G10*data.EC,data.E)+data.beta*kron(data.G11Int,data.E);

A2 = [a_ksi_2, a_theta_2;
      b_ksi_2, b_theta_2];
    end 


end




