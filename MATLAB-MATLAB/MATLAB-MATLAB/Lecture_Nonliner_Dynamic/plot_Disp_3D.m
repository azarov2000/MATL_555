function plot_Disp_3D(t,Z0,TN,tNS,m,index_disp)

tN=t/TN;    % переходим от безразмерного времени к количеству оборотов
I=find(tN>(tN(end)-tNS)); % вытаскиваем индексы последних оборотов


%% График перемещений

vect_length_beam = 0:1/m:1;

%% Формирование векторов 
%lenI = length(I);
    x_disp = [zeros(length(I),1),Z0(I(:),index_disp(1:2:length(index_disp))),zeros(length(I),1)];
    y_disp = [zeros(length(I),1),Z0(I(:),index_disp(2:2:length(index_disp))),zeros(length(I),1)];
    for j=1:length(I)
    z_disp(j,:) = vect_length_beam;
    end
    %z_disp = [zeros(length(I),1),zeros(length(I),length(1:2:11))];

%%
Ratio = [-0.030,0.030]*1;
figure('WindowState','maximized');
hold on;box on;hold on
for j=1:1:1000
    box on
    plot3(z_disp(j,:),y_disp(j,:),x_disp(j,:),'.b','MarkerSize',14)
%     FF(j) = getframe(gcf);
end
%plot3(x_disp(:,:),y_disp(:,:),z_disp(:,:),'.-r')
ff = gca; 
ff.FontName = 'Times New Roman';
ff.FontSize = 16; 

xlabel('\itz');
ylabel('\ity');
zlabel('\itx');
grid on; hold on; box on;
xlim([0 1.1]);
ylim(Ratio);
zlim(Ratio);

% vd = VideoWriter('VIDEO2');
% vd.FrameRate = 500;
% open(vd);
% for i=1:length(FF)
%    writeVideo(vd,FF(i)); 
% end
% 
% close(vd)



end 