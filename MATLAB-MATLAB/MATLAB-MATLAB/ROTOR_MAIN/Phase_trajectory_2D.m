function Phase_trajectory_2D(t,Z0,TN,tNS,disp,N)

tN=t/TN;                                                % переходим от безразмерного времени к количеству оборотов
tPoincare=(t(1):TN:t(end))';                            % вектор времени, разбитый на равные интервалы, равные времени одного оборота ротора
ZxPoincare=interp1(t,Z0(:,disp(1)),tPoincare,'spline'); % вектор перемещений по oX через каждый оборот ротора 
ZyPoincare=interp1(t,Z0(:,disp(2)),tPoincare,'spline'); % вектор перемещений по oY через каждый оборот ротора
tPoincareIndex=tPoincare/TN;                            % вектор количества оборотов       
II=find(tPoincareIndex>(tPoincareIndex(end)-tNS));     % достаём индексы точек стробирования, соотвестствующие tNS2 последним оборотам
I=find(tN>(tN(end)-tNS));                               % достаём индексы точек самой траекториии, соотв. tNS последним оборотам                                

% Построение фазовой траектории с точкам стробирования

figure('WindowState','maximized'); 
    hold on
    plot(-Z0(I,disp(2)),Z0(I,disp(1)),'LineWidth',0.5);
    plot(-ZyPoincare(II),ZxPoincare(II),'r.','MarkerSize',18);
    axis equal
    xlabel('\xi_{\it y}')
    ylabel('\xi_{\it x}','Rotation',0)
    title(['N = ',num2str(N),'; Количество последних оборотов: ',num2str(tNS)])
    ax1 = gca;
    ax1.FontName = 'Times New Roman';
    ax1.FontSize = 18;
    grid on; hold on; box on;
    xlim padded
    ylim padded





end