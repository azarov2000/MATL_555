function information_plot(dat)




figure;
axis('off');
% Вывод параметров в столбец

text(0, 0.9, ['m ', num2str(dat.m)],'FontSize', 14);
text(0, 0.8, ['N: ', num2str(dat.N)],'FontSize', 14);
text(0, 0.7, ['M_{кр}: ', num2str(dat.Mom)],'FontSize', 14);

text(0, 0.6, ['\mu_{ B}: ', num2str(dat.mu_B)],'FontSize', 14);
text(0, 0.5, ['\kappa_{ B}: ', num2str(dat.kappa_B)],'FontSize', 14);
text(0, 0.4, ['\eta_{ B}: ', num2str(dat.eta_B)],'FontSize', 14);

text(0.4, 0.9, ['\zeta_{ e}: ', num2str(dat.zeta_e)],'FontSize', 14);
text(0.4, 0.8, ['\zeta_{ V}: ', num2str(dat.zeta_V)],'FontSize', 14);


text(0.4, 0.7, ['excenty_{ Amp}: ', num2str(dat.Ampl_eps)],'FontSize', 14);
text(0.4, 0.6, ['phase_{ Amp}: ', num2str(dat.Ampl_phase)],'FontSize', 14);


end

