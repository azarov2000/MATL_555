function [ksi_x, ksi_y,ksi_x_dot,ksi_y_dot] = get_KSI(dat,Z0)

h = 1/dat.m;
%% Формирование ksi_x и ksi_y в каждом узле
flag = 1;
for j=1:2:(length(dat.A2)-2)
    for i=1:length(Z0(:,1))
        ksi_x(i,j) = Z0(i,j) + Z0(i,length(dat.A0)-1)*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_x(:,2:2:end) = [];
ksi_x = [zeros(length(ksi_x(:,1)),1),ksi_x,Z0(:,length(dat.A0)-1)];

flag = 1;
for j=2:2:(length(dat.A2)-2)
    for i=1:length(Z0(:,1))
        ksi_y(i,j) = Z0(i,j) + Z0(i,length(dat.A0))*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_y(:,1:2:end) = [];
ksi_y = [zeros(length(ksi_y(:,1)),1),ksi_y,Z0(:,length(dat.A0))];

z_vect = 0:h:1;

%% Формирование скоростей ksi_x и ksi_y в каждом узле
flag = 1;
for j=length(dat.A2)+1:2:2*(length(dat.A2))-2
    for i=1:length(Z0(:,1))
        ksi_x_dot(i,j) = Z0(i,j) + Z0(i,2*length(dat.A0)-1)*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_x_dot(:,1:length(dat.A2))= [];
ksi_x_dot(:,2:2:end) = [];
ksi_x_dot = [zeros(length(ksi_x_dot(:,1)),1),ksi_x_dot,Z0(:,2*length(dat.A0)-1)];

flag = 1;
for j=length(dat.A2)+2:2:2*(length(dat.A2))-1
    for i=1:length(Z0(:,1))
        ksi_y_dot(i,j) = Z0(i,j) + Z0(i,2*length(dat.A0))*(3*(flag*h)^2-2*(flag*h)^3);
    end
    flag = flag + 1;
end
ksi_y_dot(:,1:length(dat.A2)+1)= [];
ksi_y_dot(:,2:2:end) = [];
ksi_y_dot = [zeros(length(ksi_y_dot(:,1)),1),ksi_y_dot,Z0(:,2*length(dat.A0))];



end

