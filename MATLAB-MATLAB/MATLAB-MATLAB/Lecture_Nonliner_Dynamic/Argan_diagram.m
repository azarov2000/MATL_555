function Argan_diagram ()
  global I E zi S Z ze Ze 
  global G00 G10 G01 G11 G00Int G01Int G10Int G11Int
  global DE EC beta muR betaR h Zte 
  
  flag = 0;
  N = 0;
  step_N = 0.1;
  i = 1;
  mb_crit = 0;
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


pause(5)


    function [A0,A1,A2] = maxtix_for_Argan(N) % Нужно пересчитывать матрицы
%         global I E zi S Z ze Ze 
%         global G00 G10 G01 G11 G00Int G01Int G10Int G11Int
%         global DE EC beta muR betaR h Zte 
            aksi=kron(I,E+2*zi*N*S);
            dksi=kron(Z,E);
            ateta=kron(Z,E);
            dteta=kron(I,S-2*zi*N*E);


            bksi=kron(2*zi*I+2*ze*h*G00Int+2*Ze*G00*EC,E);
            eksi=kron(2*betaR*N*(-h*G01Int+G00*DE)-2*beta*N*G01*EC,E)+kron(2*Zte*G01*EC,S);         % поставил '-'
            bteta=kron(2*ze*h*G10Int+2*Ze*G10*EC,E);
            eteta=kron(2*betaR*N*(-h*G11Int+G10*DE)-2*beta*N*G11*EC,E)+kron(2*zi*I+2*Zte*G11*EC,S); % поставил '-'


            cksi=kron(muR*h*G00Int+G00*EC,E);
            fksi=kron(betaR*(h*G01Int-G00*DE)+beta*G01*EC,S);
            cteta=kron(muR*h*G10Int+G10*EC,E);
            fteta=kron(betaR*(h*G11Int-G10*DE)+beta*G11*EC,S);


            A0=[aksi dksi;ateta dteta];

            A1=[bksi eksi;bteta eteta];

            A2=[cksi fksi; cteta fteta];
    end 


end




