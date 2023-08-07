function Spectral_analysis_of_displacements(t,Z0,tNS,TN,disp)


tN=t/TN;                    % переходим от безразмерного времени к количеству оборотов
I=find(tN>(tN(end)-tNS));   % вытаскиваем индексы последних оборотов

tend=(tN(I)*TN);            % время последних оборотов

Zx=(Z0(I,disp(1)))';
Zy=Z0(I,disp(2))';
ttInt=linspace(tend(1),tend(end),length(tend));
ZZx=interp1(tend,Zx,ttInt);
ZZy=interp1(tend,Zy,ttInt);
PowerSignal=0.999;         % процент мощности возращаемого сигнала
[fx,Xspectr]=spectrum_Fig_Power(ttInt,ZZx,0,1,PowerSignal);
[fy,Yspectr]=spectrum_Fig_Power(ttInt,ZZy,0,1,PowerSignal);
%%
figure('WindowState','maximized');
     subplot(2,2,1)
     grid on; hold on; box on;
     plot(ttInt,ZZx)
     xlabel('\tau')
     ylabel('\xi_{\itx}')
     title(['Рассматривается ',num2str(tNS),' последних оборотов'])
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;
     xlim padded
     ylim padded
     
     subplot(2,2,2)
     h1=stem(fx,Xspectr);
     grid on; hold on; box on;
     xline(1/TN,'--r',[{'Частота вращения ротора'};'(',num2str(1/TN),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('\it f')
     ylabel('Амплитуда')
     title('Спектральный анализ \xi_{\it x}')
     grid on;
     set(get(h1,'BaseLine'),'LineStyle','-');
     set(h1,'MarkerFaceColor','red');
     xlim padded
     ylim padded
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;
     grid on;
     
     subplot(2,2,3)
     grid on; hold on; box on;
     plot(ttInt,ZZy)
     xlabel('\tau')
     ylabel('\xi_{\ity}')
     title(['Рассматривается ',num2str(tNS),' последних оборотов'])
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;
     xlim padded
     ylim padded
     
     subplot(2,2,4)
     grid on; hold on; box on;
     h2=stem(fy,Yspectr);
     xline(1/TN,'--r',[{'Частота вращения ротора'};'(',num2str(1/TN),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('\it f')
     ylabel('Амплитуда')
     title('Спектральный анализ \xi_{\it y}')
     set(get(h2,'BaseLine'),'LineStyle','-');
     set(h2,'MarkerFaceColor','red');
     xlim padded
     ylim padded
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;

end