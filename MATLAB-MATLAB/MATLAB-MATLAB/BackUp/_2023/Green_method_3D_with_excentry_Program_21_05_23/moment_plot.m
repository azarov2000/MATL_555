function moment_plot(N_end,t_m,Mi,dat,sol)

index = find(t_m>sol.TN*N_end);
n_rotation = t_m(index)/sol.TN;

figure;
box on; grid on; hold on
plot(n_rotation,Mi(index));
xlabel('Кол-во оборотов')
ylabel('M_{torsion}')
ff = gca;
ff.FontSize = 16;

title(['M = ',num2str(dat.Mom),'N = ',num2str(dat.N)]);




end

