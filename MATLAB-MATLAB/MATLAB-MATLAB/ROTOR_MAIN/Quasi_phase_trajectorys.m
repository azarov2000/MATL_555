function Quasi_phase_trajectorys(t,Z0,TN,tNS,disp)

    tEta=t/TN;                                                   % Переходим к количеству оборотов
    NumerOfRev=find(tEta>(tEta(end)-tNS));                       % Вытаскиваем индексы количества оборотов
    [~,MaxKSIx]=ext(tEta(NumerOfRev),Z0(NumerOfRev,disp(1)));
    [~,MaxKSIy]=ext(tEta(NumerOfRev),Z0(NumerOfRev,disp(2)));

    vect_xx = MaxKSIx(1:(length(MaxKSIx)-1));
    vect_xy = MaxKSIx(2:length(MaxKSIx));
    
    vect_yx = MaxKSIy(1:(length(MaxKSIy)-1));
    vect_yy = MaxKSIy(2:length(MaxKSIy));
    
   
    figure;
        subplot(121)
        grid on; hold on; box on;
        p = plot(vect_xx(:),vect_xy(:),'.');
        p.MarkerSize = 20;
        p.MarkerFaceColor = 'k';
        x = [-0.5 ,  0.5];
        y = [-0.5 ,  0.5];
        pl = line(x,y);
        pl.Color = 'b';
        pl.LineStyle = '--';
        title(['Квазифазовая траектория за ',num2str(tNS),' последних оборотов'] )
        ff = gca; 
        ff.FontName = 'Times New Roman';
        ff.FontSize = 20;
        xlabel('\xi_{\it x}^{extr}(n)');
        ylabel('\xi_{\it x}^{extr}(n+1)');
        daspect([1 1 1])
        
        subplot(122)
        grid on; hold on; box on;
        p = plot(vect_yx(:),vect_yy(:),'.');
        p.MarkerSize = 20;
        p.MarkerFaceColor = 'k';
        x = [-0.5 ,  0.5];
        y = [-0.5 ,  0.5];
        pl = line(x,y);
        pl.Color = 'b';
        pl.LineStyle = '--';
        title(['Квазифазовая траектория за ',num2str(tNS),' последних оборотов'] )
        ff = gca; 
        ff.FontName = 'Times New Roman';
        ff.FontSize = 20;
        xlabel('\xi_{\it y}^{extr}(n)');
        ylabel('\xi_{\it y}^{extr}(n+1)');
        daspect([1 1 1])


end