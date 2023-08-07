function plot_Disp_3D(t,Z0,TN,tNS,m,index_disp,jC,nameS)

tN=t/TN;    % переходим от безразмерного времени к количеству оборотов
I=find(tN>(tN(end)-tNS)); % вытаскиваем индексы последних оборотов

% вектор координат узлов коллокации
vect_length_beam = 0:1/m:1;

switch nameS
    case 1
%% Формирование векторов 
    x_disp = [zeros(length(I),1),Z0(I(:),index_disp(1:2:length(index_disp))),zeros(length(I),1)];
    y_disp = [zeros(length(I),1),Z0(I(:),index_disp(2:2:length(index_disp))),zeros(length(I),1)];
    for j=1:length(I)
    z_disp(j,:) = vect_length_beam;
    end
    case 2
%% Формирование векторов 
    x_disp = [zeros(length(I),1),Z0(I(:),index_disp(1:2:length(index_disp)))];
    y_disp = [zeros(length(I),1),Z0(I(:),index_disp(2:2:length(index_disp)))];
    for j=1:length(I)
    z_disp(j,:) = vect_length_beam;
    end
end 
   
%% График перемещений
Ratio = [-0.030,0.030]*1;
fig1 = figure('WindowState','maximized');
axes1 = axes('Parent',fig1);
    hold on;box on;hold on

    for j=1:length(I) % цикл по всем 
        box on
        plot3(z_disp(j,[1:1:jC,jC+2:end]),y_disp(j,[1:1:jC,jC+2:end]),x_disp(j,[1:1:jC,jC+2:end]),'.b','MarkerSize',14)
        plot3(z_disp(j,jC+1),y_disp(j,jC+1),x_disp(j,jC+1),'.r','MarkerSize',14)
    end
    for j=1:length(z_disp(1,:))
        p = plot3(z_disp(1,j),0,0,'ok','MarkerSize',10);
        p.MarkerFaceColor = [1 0.5 0];
    end
    
     pl = line([1,0],[0,0],[0,0],'LineWidth',2,'LineStyle','--');
     pl.Color = 'k';
     title(['Количество последних оборотов: ',num2str(tNS)])

    view(axes1,[24.237499703950355,33.634394598615501]);
    ff = gca; 
    ff.FontName = 'Times New Roman';
    ff.FontSize = 20; 

    xlabel('\zeta');
    ylabel('\it y');
    zlabel('\it x','Rotation',0);
    grid on; hold on; box on;
    xlim([0 1.1]);
    ylim(Ratio);
    zlim(Ratio);
end 