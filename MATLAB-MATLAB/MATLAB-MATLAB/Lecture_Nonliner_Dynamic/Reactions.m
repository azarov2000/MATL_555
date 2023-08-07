function Reactions(t,Z0,TN,tNS,NumbOp,disp,disp_vel,kappa_s,eta_s,D_s,f_s,N,eD)


tN=t/TN;                  % переходим от безразмерного времени к количеству оборотов
I=find(tN>(tN(end)-tNS)); % вытаскиваем индексы последних оборотов

% Запись сил в упругой опоре
ltt=length(tN); PItog=zeros(2,ltt); flag=1;
 
for alfa1=NumbOp
     for k=1:ltt 
        alfa=(alfa1*pi)/180; % переход к радианам
        % Связь перемещений
        uRxs=Z0(k,disp(1))*cos(alfa)+Z0(k,disp(2))*sin(alfa);
        uRys=-Z0(k,disp(1))*sin(alfa)+Z0(k,disp(2))*cos(alfa);
        % Связь скоростей
        duRxs=Z0(k,disp_vel(1))*cos(alfa)+Z0(k,disp_vel(2))*sin(alfa);
        duRys=-Z0(k,disp_vel(1))*sin(alfa)+Z0(k,disp_vel(2))*cos(alfa);
        % Связь сил
        Sn(flag,k)=(kappa_s*(uRxs-eta_s)+D_s*duRxs)*((uRxs-eta_s)>=0)*((kappa_s*(uRxs-eta_s)+D_s*duRxs)>=0);
        St(flag,k)=f_s*((N*eD/2+duRys)/abs(N*eD/2+duRys))*Sn(flag,k);
        MatrRot=[cos(alfa), -sin(alfa); sin(alfa),cos(alfa)];
        PItog(:,k)=MatrRot*[-Sn(flag,k);-St(flag,k)];
     end
     flag=flag+1;
end

 
 %% Графики сил в упругих опорах  опоре

lenOp=length(NumbOp); % Автоматическая запись легенд
Legend{1}=zeros(1,lenOp);
for j=1:1:lenOp
    Legend{j}=['Опора №', num2str(j)];
end

% Графики сил от действия упругих опор полном промежутке 
figure;
    subplot(221);hold on;box on;grid on
       for j=1:1:lenOp
           plot(tN,Sn(j,:))
       end
       ff = gca; 
       ff.FontName = 'Times New Roman';
       ff.FontSize = 20;
       legend(Legend)
       title('Нормальные реакции в опорах')
       xlabel('Количество оборотов')
       ylabel('S_ { n\alpha_{j}}')
    
    subplot(222);hold on;box on;grid on
       for j=1:1:lenOp
           plot(tN,St(j,:))
       end
       ff = gca; 
       ff.FontName = 'Times New Roman';
       ff.FontSize = 20;
        legend(Legend)
        title('Тангенциальные реакции в опорах')
        xlabel('{\itt} / (2 \pi / N )')
        xlabel('Количество оборотов')
        ylabel('S_ { t\alpha_{j}}')
        
    subplot(223);hold on;box on;grid on
       for j=1:1:lenOp
           plot(tN(I),Sn(j,I))
       end
       ff = gca; 
       ff.FontName = 'Times New Roman';
       ff.FontSize = 20;
       legend(Legend)
       title('Нормальные реакции в опорах')
       xlabel('Количество оборотов')
       ylabel('S_ { n\alpha_{j}}')
       
    subplot(224);hold on;box on;grid on
       for j=1:1:lenOp
           plot(tN(I),St(j,I))
       end
       ff = gca; 
       ff.FontName = 'Times New Roman';
       ff.FontSize = 20;
        legend(Legend)
        title('Тангенциальные реакции в опорах')
        xlabel('{\itt} / (2 \pi / N )')
        xlabel('Количество оборотов')
        ylabel('S_ { t\alpha_{j}}')
end