function Precession_coefficient(t,Z0,TN,tNS,disp,N)

tN=t/TN;    % переходим от безразмерного времени к количеству оборотов
I=find(tN>(tN(end)-tNS)); % вытаскиваем индексы последних оборотов

timePI = TN*tN(I); % время, соответствующее последним оборотам

%% Измнение центра траектории
Px=polyfit(timePI,Z0(I,disp(1)),0);
Py=polyfit(timePI,Z0(I,disp(2)),0);
TrandX=polyval(Px,timePI);    % Получаем записываем ординаты этого полинома
TrandY=polyval(Py,timePI);    % Получаем записываем ординаты этого полинома
 
Z0(I,disp(1))=Z0(I,disp(1))-TrandX; % по сути мы убираем постоянную составляющую сигнала
Z0(I,disp(2))=Z0(I,disp(2))-TrandY;

%%




[tPI,PI] = PrecessionCoeff(timePI,Z0(I,disp(1)),Z0(I,disp(2)),N);
PowerSignal = 0.999;
[fPi,RealizationPi] = spectrum_Fig_Power(tPI,PI,0,PowerSignal);
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
         ylim([-1,1.5])
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
         xlim([0,6])
         grid on;
     

     

end