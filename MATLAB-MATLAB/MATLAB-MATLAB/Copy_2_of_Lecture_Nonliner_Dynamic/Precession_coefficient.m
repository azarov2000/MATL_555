function Precession_coefficient(t,Z0,TN,tNS,disp,N)

tN=t/TN;    % переходим от безразмерного времени к количеству оборотов
I=find(tN>(tN(end)-tNS)); % вытаскиваем индексы последних оборотов

timePI = TN*tN(I); % время, соответствующее последним оборотам
[tPI,PI] = PrecessionCoeff(timePI,Z0(I,disp(1)),Z0(I,disp(2)),N);
PowerSignal = 0.999;
[fPi,RealizationPi] = spectrum_Fig_Power(tPI,PI,0,1,PowerSignal);
Average = mean(PI);

figure('WindowState','maximized');
      subplot(2,1,1)
         plot(tPI,PI,'LineWidth',1)
         xlabel('\tau','FontName','Times New Roman','FontSize',16)
         ylabel('\Lambda','FontName','Times New Roman','FontSize',16)
         text((max(tPI)+min(tPI))/2,(0+1)/2,['Avg(\Lambda) = ',num2str(Average)],...
             'FontName','Times New Roman','FontSize',20,'BackGround','white','EdgeColor','black')
         title(['Количество последних оборотов: ',num2str(tNS)])
         ff = gca; 
         ff.FontName = 'Times New Roman';
         ff.FontSize = 20;
         xlim padded
         ylim([0,1.5])
         grid on;
     
     subplot(2,1,2)
         h1=stem(fPi,RealizationPi);
         xline(1/TN,'--r',[{'Частота вращения ротора'};'(',num2str(1/TN),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
         xlabel('Безразмерная частота','FontName','Times New Roman','FontSize',16)
         ylabel('Амплитуда','FontName','Times New Roman','FontSize',16)
         title(['Спектр \Lambda, мощность: ',num2str(PowerSignal*100),'%'])
         grid on;
         set(get(h1,'BaseLine'),'LineStyle','-');
         set(h1,'MarkerFaceColor','red');
         ff = gca; 
         ff.FontName = 'Times New Roman';
         ff.FontSize = 20;
         grid on;
     
     function [tPr,CoeffPr] = PrecessionCoeff(t_init, x_init,y_init,N)
           length_t=length(t_init);
           tInit=linspace(t_init(1),t_init(end),length_t);
           xInit=interp1(t_init,x_init,tInit);
           yInit=interp1(t_init,y_init,tInit);

            h=(tInit(end)-tInit(1))/(length_t-1);
            tPr=tInit(2:(end-1));
              X=xInit(2:(end-1));
              Y=yInit(2:(end-1));
              D=X.^2+Y.^2;
                DX=xInit(3:end)-xInit(1:(end-2));
                DY=yInit(3:end)-yInit(1:(end-2));           
                         fr=(DY.*X-DX.*Y)./D;
                    CoeffPr=fr/(2*h*N);

     end
     

end