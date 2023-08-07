function plot_3D(N_end,t,ksi_x,ksi_y,dat,sol)

index = find(t>=sol.TN*N_end);
h = 1/dat.m;
z_vect = 0:h:1;
for i=1:length(z_vect)
    Z_3D{i} = kron(ones(length(ksi_x),1),z_vect(i));
end
figure;
hold on; box on; grid on;
view([45.9000000862685 13.7999998741406])
for i=1:length(z_vect)
    plot3(Z_3D{i}(index),ksi_x(index,i),ksi_y(index,i),'.-k')
end
xlabel('\zeta')
ylabel('\xi_{ x}')
zlabel('\xi_{ y}','Rotation',0)
ff = gca;
ff.FontSize = 18;
title(['M = ',num2str(dat.Mom),', ','N = ',num2str(dat.N)]);


end

