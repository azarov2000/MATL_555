function disp_param(d)

current_time = clock;

fig = figure;
set(fig, 'Name', 'Параметры системы');
text(0.3,1,'Параметры системы','FontSize', 12, 'FontWeight', 'bold')

text(0.33,0.95,['(Время запуска: ',num2str(current_time(4)),':',num2str(current_time(5)),')'],'FontSize', 10)

text(0.01,0.8,['m = ',num2str(d.m)],'FontSize', 14)
text(0.01,0.7,['\eta_e = ',num2str(d.zeta_e)],'FontSize', 14)
text(0.01,0.6,['\eta_V = ',num2str(d.zeta_V)],'FontSize', 14)
text(0.01,0.5,['\eta_{VV} = ',num2str(d.zeta_VV)],'FontSize', 14)
text(0.01,0.4,['\alpha = ',num2str(d.alpha)],'FontSize', 14)
text(0.01,0.3,['\Gamma_{0} = ',num2str(d.Gamma_0)],'FontSize', 14)

text(0.4,0.8,'Расчётные параметры','FontSize', 14,'FontWeight','bold')
text(0.4,0.7,['\Omega_{crit} = ',num2str(d.Ncritical)],'FontSize', 14)
text(0.4,0.6,['\nu_{1} = ',num2str(d.natural_freq(1))],'FontSize', 14)
text(0.4,0.5,['\nu_{2} = ',num2str(d.natural_freq(2))],'FontSize', 14)
text(0.4,0.4,['\nu_{3} = ',num2str(d.natural_freq(3))],'FontSize', 14)
text(0.4,0.3,['\nu_{biff} = ',num2str(d.freq_biff)],'FontSize', 14)

axis off

end

