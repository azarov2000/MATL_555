function Spectral_analysis(d)

index = find(d.t >= d.tNS);

t_Int = linspace(d.t(1),d.t(end),length(d.t));
ksi_x_Int = interp1(d.t,d.ksi_x,t_Int);
ksi_y_Int = interp1(d.t,d.ksi_y,t_Int);
[fx,Xspectr]=spectrum_Fig_Power(t_Int(index),ksi_x_Int(index),0,d.Power);
[fy,Yspectr]=spectrum_Fig_Power(t_Int(index),ksi_y_Int(index),0,d.Power);
%%
figure('WindowState','maximized');
     subplot(2,2,1)
     grid on; hold on; box on;
     plot(t_Int,ksi_x_Int)
     xlabel('\tau')
     ylabel('\xi_{\it x}')
     title(['Рассматривается ',num2str(d.tNS),' последних секунд'])
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;
     xlim padded
     ylim padded
     
     subplot(2,2,2)
     h1=stem(fx,Xspectr);
     grid on; hold on; box on;
     xline(d.N/(2*pi),'--r',[{'Частота вращения ротора'};'(',num2str(d.N/(2*pi)),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
%      xline(d.nat_freq(1)/(2*pi),'--r',[{'Первая СЧ'};'(',d.nat_freq(1)/(2*pi),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
%      xline(d.nat_freq(2)/(2*pi),'--r',[{'Вторая СЧ'};'(',d.nat_freq(2)/(2*pi),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
%      xline(d.nat_freq(3)/(2*pi),'--r',[{'Третья СЧ'};'(',d.nat_freq(3)/(2*pi),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('Безразмерная частота')
     ylabel('Амплитуда')
     title(['Спектр \xi_{\it x}; мощность ',num2str(d.Power*100),'%'])
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
     plot(t_Int,ksi_y_Int)
     xlabel('\tau')
     ylabel('\xi_{\it y}')
     title(['Рассматривается ',num2str(d.tNS),' последних секунд'])
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;
     xlim padded
     ylim padded
     
     subplot(2,2,4)
     grid on; hold on; box on;
     h2=stem(fy,Yspectr);
     xline(d.N/(2*pi),'--r',[{'Частота вращения ротора'};'(',num2str(d.N/(2*pi)),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
%      xline(d.nat_freq(1)/(2*pi),'--r',[{'Первая СЧ'};'(',d.nat_freq(1)/(2*pi),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
%      xline(d.nat_freq(2)/(2*pi),'--r',[{'Вторая СЧ'};'(',d.nat_freq(2)/(2*pi),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
%      xline(d.nat_freq(3)/(2*pi),'--r',[{'Третья СЧ'};'(',d.nat_freq(3)/(2*pi),')'],'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('Безразмерная частота')
     ylabel('Амплитуда')
     title(['Спектр \xi_{\it y}; мощность ',num2str(d.Power*100),'%'])
     set(get(h2,'BaseLine'),'LineStyle','-');
     set(h2,'MarkerFaceColor','red');
     xlim padded
     ylim padded
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;

end