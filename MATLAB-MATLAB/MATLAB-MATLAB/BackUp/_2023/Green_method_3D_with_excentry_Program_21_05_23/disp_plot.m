function disp_plot(N_end,t,ksi_x,ksi_y,dat,sol)

index = find(t>=sol.TN*N_end);
n_rotation = t(index)/sol.TN;

figure;
subplot(2,1,1)
hold on; box on; grid on;
plot(n_rotation,ksi_x(index, end))
xlabel('Кол-во оборотов')
ylabel('\xi_{ x}','Rotation', 0)
ff = gca; ff.FontSize = 18;
subplot(2,1,2)
hold on; box on; grid on;
plot(n_rotation,ksi_y(index, end))
xlabel('Кол-во оборотов')
ylabel('\xi_{ y}','Rotation', 0)
ff = gca; ff.FontSize = 18;

%suptitle(['M = ',num2str(dat.Mom),'N = ',num2str(dat.N)]);


end

