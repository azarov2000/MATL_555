function Spectral_analysis(t,Z0,tNS,TN,disp)


tN=t/TN;                    % переходим от безразмерного времени к количеству оборотов
I=find(tN>(tN(end)-tNS));   % вытаскиваем индексы последних оборотов

tend=(tN(I)*TN);            % время последних оборотов

Zx=(Z0(I,disp(1)));
Zy=Z0(I,disp(2));
ttInt=linspace(tend(1),tend(end),length(tend));
ZZx=interp1(tend,Zx,ttInt);
ZZy=interp1(tend,Zy,ttInt);
PowerSignal=0.999;         % процент мощности возращаемого сигнала
[fx,Xspectr]=spectrum_Fig_Power(ttInt,ZZx,0,1,PowerSignal);

% Поиск эксремумов 

[text,yext] = ext(ttInt,ZZx);




%%
figure('WindowState','maximized');
     subplot(1,2,1)
     grid on; hold on; box on;
     plot(ttInt,ZZx)
     plot(text,yext,'.r','MarkerSize',20)
     xlabel('\tau')
     ylabel('\xi_{\itx}')
     title(['Рассматривается ',num2str(tNS),' последних оборотов'])
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;
     xlim padded
     ylim padded
     ll = length(yext);
     subplot(1,2,2)
     xx = yext(1:2:ll-1);
     yy = yext(2:2:ll);
     for i = 1:length(xx)
        plot(xx(i),yy(i),'.r','MarkerSize',20);
        hold on;
     end
         
     xlabel('\xi^{\it extr}_{\it x}(n)')
     ylabel('\xi^{\it extr}_{\it x}(n+1)')
     xlim padded
     ylim padded
     ff = gca;
     ff.FontName = 'Times New Roman';
     ff.FontSize = 18;
     grid on;
     

end