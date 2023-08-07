function Phase_trajectory_2D(t,Z0,TN,nT,tNS,tNS2,disp,N)

tN=t/TN;                                                % переходим от безразмерного времени к количеству оборотов
tPoincare=(t(1):TN:t(end))';                            % вектор времени, разбитый на равные интервалы, равные времени одного оборота ротора
ZxPoincare=interp1(t,Z0(:,disp(1)),tPoincare,'spline'); % вектор перемещений по oX через каждый оборот ротора 
ZyPoincare=interp1(t,Z0(:,disp(2)),tPoincare,'spline'); % вектор перемещений по oY через каждый оборот ротора
tPoincareIndex=tPoincare/TN;                            % вектор количества оборотов       
II=find(tPoincareIndex>(tPoincareIndex(end)-tNS2));     % достаём индексы точек стробирования, соотвестствующие tNS2 последним оборотам
I=find(tN>(tN(end)-tNS));                               % достаём индексы точек самой траекториии, соотв. tNS последним оборотам                                

% Построение фазовой траектории с точкам стробирования

figure('WindowState','maximized'); 
    hold on
    plot(-Z0(I,disp(2)),Z0(I,disp(1)),'LineWidth',0.5);
    lenPoin=length(II);
    for j=1:1:lenPoin
        if j<=3
            kColor='o';
            kSize=14;
            p=plot(-ZyPoincare(II(j)),ZxPoincare(II(j)),kColor,'MarkerSize',kSize);
            p.MarkerFaceColor = 'k';
            p.MarkerSize = 8;
        else
            kColor='r.';
            kSize=18;
            plot(-ZyPoincare(II(j)),ZxPoincare(II(j)),kColor,'MarkerSize',kSize);
        end
    end
    axis equal
    xlabel('\xi_{\it y}')
    ylabel('\xi_{\it x}','Rotation',0)
    title([num2str(nT-tNS2),' < nT < ',num2str(nT),'; ','N = ',num2str(N)])
    ax1 = gca;
    ax1.FontName = 'Times New Roman';
    ax1.FontSize = 18;
    grid on; hold on; box on;
    xlim padded
    ylim padded





end