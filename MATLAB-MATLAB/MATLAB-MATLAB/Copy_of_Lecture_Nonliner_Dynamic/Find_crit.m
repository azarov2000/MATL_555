function [Ncritical,MaxOm,MinOm] = Find_crit (N,step_N,eta_i)
  global I E S Z eta_e eta_Dksi 
  global G00 G10 G01 G11 G00Int G01Int G10Int G11Int
  global DE EC beta muR betaR h eta_Dte 
  
  flag = 0;    % флаг отсчёта 10 красных точек
  i = 1;       % счётчик цикла while
  mb_crit = 0; % начальное приближение для определения критической скорости
  while (flag<10)
      [A0,A1,A2] = maxtix_for_Argan(N,eta_i);
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
    maxOm(j) = max(imag(EigenValues{j}));
    minOm(j) = min(imag(EigenValues{j}));
end
ZEROFUN=@(x) interp1(N_vect(:),ReMax(:),x,'spline');
Ncritical=fzero(ZEROFUN,mb_crit);
MaxOm = interp1(N_vect(:),maxOm(:),Ncritical,'spline');
MinOm = interp1(N_vect(:),minOm(:),Ncritical,'spline');
%disp(['Критическая скорость равна ',num2str(Ncritical)]);


% Функция возвращает матрицы коэффициентов при заданном N
    function [A0,A1,A2] = maxtix_for_Argan(N,eta_i) 
            aksi=kron(I,E+2*eta_i*N*S);
            dksi=kron(Z,E);
            ateta=kron(Z,E);
            dteta=kron(I,S-2*eta_i*N*E);

            bksi=kron(2*eta_i*I+2*eta_e*h*G00Int+2*eta_Dksi*G00*EC,E);
            eksi=kron(2*betaR*N*(-h*G01Int+G00*DE)-2*beta*N*G01*EC,E)+kron(2*eta_Dte*G01*EC,S);
            bteta=kron(2*eta_e*h*G10Int+2*eta_Dksi*G10*EC,E);
            eteta=kron(2*betaR*N*(-h*G11Int+G10*DE)-2*beta*N*G11*EC,E)+kron(2*eta_i*I+2*eta_Dte*G11*EC,S);


            cksi=kron(muR*h*G00Int+G00*EC,E);
            fksi=kron(betaR*(h*G01Int-G00*DE)+beta*G01*EC,S);
            cteta=kron(muR*h*G10Int+G10*EC,E);
            fteta=kron(betaR*(h*G11Int-G10*DE)+beta*G11*EC,S);


            A0=[aksi dksi;ateta dteta];

            A1=[bksi eksi;bteta eteta];

            A2=[cksi fksi; cteta fteta];
    end 


end