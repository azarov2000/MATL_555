%%
close all
clc
clear
%% Дано
m = 6;          % Количество участков
mnew = m-2;
tNS=300;        % Количество оборотов для фазовой трактории
tNS2=300;       % Количество оборотов для стробоскопического отображения
%ii=6;           % Для m=6; Zс=0.7; (ii=0), если диск расположен посередине
ii=2;
eta_s = (10*10^-3)/0.7;  % коэффициент зазора

zDisk=0.7;        % Кооридната расположения диска
h=1/m;          % Граничные случаи задание постановки диска
s=round(zDisk/h);
if s==0
    jC=1;
elseif s==m
    jC=m;
else
    jC=s;
end
% Вектор начальных условий
x01=zeros((mnew+1),1); x01(jC,1)=1; x0=kron(x01,ones(2,1));
x0=0.1*0.01*[x0;zeros(2*(mnew+1),1);zeros(4*(mnew+1),1)];

%% Вектора 
vect_N = [14.5:0.25:20,20.5:0.5:45];
le_vect_N = length(vect_N);
%% Решение системы дифференциальных уравнений
x01=zeros((mnew+1),1); x01(jC,1)=1; x0=kron(x01,ones(2,1));
x0=0*[x0;zeros(2*(mnew+1),1);zeros(4*(mnew+1),1)];
% vect_N = [14.5:0.25:20,20.5:0.5:50];
vect_N = 11:0.25:13.5;
le_vect_N = length(vect_N);


NUMREV = 300;
tic;
for j=1:le_vect_N
    j
    [t,Z0,t_for_spectrum,Coeff_for_spectrum,Average] = FINDSOLVE(vect_N(j),x0,NUMREV);
    x0 = Z0(end,:);
    
    tayPI{j} = t_for_spectrum;
    PI{j} = Coeff_for_spectrum;
    AveragePi(j) = Average
    % Сохранение реализаций
    if (j==1||j==2||j==3||j==4||j==5||j==6||j==7||j==8||j==9||j==10||j==11)
        ZZ{j} = Z0; time{j} = t;
    end
    clear t Z0 tPr CoeffPr Average fPi RealizationPi Trand
end
toc;

le_Last_realization = length(tayPI{le_vect_N});
for j=74:1:le_vect_N-1
    tayPI{j} = tayPI{j}(1:le_Last_realization);
    PI{j} = PI{j}(1:le_Last_realization);
end
%%
PowerSignal = 1;
for j=1:1:le_vect_N
    
    [fPi,RealizationPi,Trand] = spectrum_Fig_Power(tayPI{j},PI{j},0,1,PowerSignal);
    fPI{j} = fPi;
    AmplitudePI{j} = RealizationPi;
    TrandD{j} = Trand;
end 
     
%% Просмотр реализаций
IND = 1;     % Индекс реализации
TN = 2*pi/vect_N(IND);    % время одного периода
nT = 2000;      % число периодов - нужно много периодов для получения установившегося решения
tN = time{IND}/TN;
tlim = TN*nT;
I = find(tN>(tN(end)-tNS)); % Индексы вектора последних оборотов
     
%% {N,Пmean}
figure('WindowState','maximized');
box on; grid on; hold on;
for j=1:1:le_vect_N
%     plot(vect_N(j)/(2*pi),AveragePi(j),'r.','MarkerSize',18)
    plot(vect_N(j),AveragePi(j),'r.','MarkerSize',18)
end
plot(vect_N(:),AveragePi(:),'-r')
% xline(18.2624/(2*pi),'--',{'N_{critical}'},'FontName','Times New Roman','FontSize',16,'LabelVerticalAlignment','middle')
xline(13.9185,'--',{'N_{critical}'},'FontName','Times New Roman','FontSize',16,'LabelVerticalAlignment','middle')
str = '$$ \overline{\Lambda} $$';
xlabel('N')
ylabel(str,'interpreter','latex')
title([num2str(nT-tNS),' < nT < ',num2str(nT)])
    ff = gca; 
    ff.FontName = 'Times New Roman';
    ff.FontSize = 20;
ylim([-1.5,1.5])

     
     %% {f,N,A_Пmean}
    figure1 = figure('WindowState','maximized');
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
     box on; grid on; hold on;
     for j=1:1:le_vect_N
         le_fPI = length(fPI{j});
         NN = ones(1,le_fPI)*vect_N(j);
         plot3(NN/(2*pi),fPI{j},AmplitudePI{j},'LineWidth',1.5,'Color','b')
     end
     %plot3([0 0 4 4 0], (14.2624/(2*pi))*[1 1 1 1 1], [0 1 1 0 0] )
     %fill3([1.7 1.7 3.25 3.25 1.7], (14.2624/(2*pi))*[1 1 1 1 1], [0 0.5 0.5 0 0],'yellow','FaceAlpha',0.3)
     title(['Cпектры \Lambda','; ',num2str(nT-tNS),' < nT < ',num2str(nT)],'FontName','Times New Roman','FontSize',16)
     text(0,0,0,'20')
     xlabel('N, 1/s')
     ylabel('f, 1/s')
     zlabel('Amplitude')
     ylim([0,20])
     len_TrendD = length(TrandD);
     for j=1:len_TrendD
         texx{j} = ['\rm При N = ',num2str(vect_N(j)),':  ', num2str(mean(TrandD{j}))];
     end
%      text(2,15,0.5,texx)
    text('Parent',axes1,'BackgroundColor','w','FontName','Times New Roman','FontSize',10,'String',...
        ['\bf Нулевая гармоника ','\bf спектрального разложения: ',' ',texx],'Position',[2 15 0.5])
    set(axes1,'FontName','Times New Roman','FontSize',18);
    
    ff = gca; 
    ff.FontName = 'Times New Roman'
    ff.FontSize = 20
 %% Графики перемещений и углов поворота
figure2 = figure('WindowState','maximized');
axes2 = axes('Parent',figure2);
subplot(211);hold on;box on
plot(tN,ZZ{IND}(:,mnew+ii+1),tN,ZZ{IND}(:,mnew+ii+2))
    
yline(-eta_s,'--k',{'support № 2'},'FontName','Times New Roman','FontSize',14);
yline(eta_s*0.5,'--k',{'supports № 1,3'},'FontName','Times New Roman','FontSize',14);
legend('\xi_{ x}','\xi_{ y}','FontName','Times New Roman','FontSize',20)
xlabel(['Displacement, ','{\itt} / (2 \pi / N )'],'FontName','Times New Roman','FontSize',16)
    
grid on;
    
%     title(['\epsilon = ',num2str(vect_disbalance(ind)),'; ',...
%         '\kappa = ' ,num2str(vect_stiff(IND))])
    
subplot(212);hold on;box on
plot(tN,ZZ{IND}(:,2*(mnew+1)+(mnew)+ii+1),tN,ZZ{IND}(:,2*(mnew+1)+(mnew)+ii+2))
    
legend('\vartheta_{ x}','\vartheta_{ y}','FontName','Times New Roman','FontSize',16)
xlabel(['Angular movements, ','{\itt} / (2 \pi / N )'],'FontName','Times New Roman','FontSize',16)
grid on;

%%
figure;
hold on; grid on; 
plot(ZZ{IND}(I,2*(mnew+1)+(mnew)+ii+2),ZZ{IND}(I,mnew+ii+1))



    
    %% Фазовая траектория
% Стробоскопическое отображение
tPoincare=(0:TN:tlim)';                                               % вектор времени, разбитый на равные интервалы, равные времени одного оборота
ZxPoincare=interp1(time{IND},ZZ{IND}(:,mnew+ii+1),tPoincare,'spline');   % Вектор перемещений по oX через каждый оборот 
ZyPoincare=interp1(time{IND},ZZ{IND}(:,mnew+ii+2),tPoincare,'spline'); % Вектор перемещений по oY через каждый оборот
tPoincareIndex=tPoincare/TN;                                          % Вектор количеств оборотов                   
II=find(tPoincareIndex>(tPoincareIndex(end)-tNS2));                   % Индексы вектора последних оборотов

%% Фазовая траектория
figure3 = figure('WindowState','maximized');
 axes3 = axes('Parent',figure3);
    hold on;box on
    plot(-ZZ{IND}(I,mnew+ii+2),ZZ{IND}(I,mnew+ii+1),'LineWidth',2);
    lenPoin=length(II);
    for j=1:1:lenPoin
        if j<=3
            kColor='o';
            kSize=14;
            p=plot(-ZyPoincare(II(j)),ZxPoincare(II(j)),kColor,'MarkerSize',kSize);
            p.MarkerFaceColor = 'k'
            p.MarkerSize = 8;
        else
            kColor='r.';
            kSize=18;
            plot(-ZyPoincare(II(j)),ZxPoincare(II(j)),kColor,'MarkerSize',kSize);
        end
    end
    
    axis equal
    xlabel('\xi_{ y}')
    ylabel('\xi_{ x}','Rotation',0)
    
%     title(['\epsilon = ',num2str(vect_disbalance(ind)),'; ',...
%         '\kappa = ' ,num2str(vect_stiff(IND))])
    grid on;
%% 
%% Спектральный анализ для одной из реализаций
tt=(tN(I)*TN);
Zx=(ZZ{IND}(I,mnew+ii+1))';
Zy=ZZ{IND}(I,mnew+ii+2)';
ttInt=linspace(tt(1),tt(end),length(tt));
ZZx=interp1(tt,Zx,ttInt);
ZZy=interp1(tt,Zy,ttInt);
PowerSignal=0.999;
[fx,Xspectr]=spectrum_Fig_Power(ttInt,ZZx,0,1,PowerSignal);
[fy,Yspectr]=spectrum_Fig_Power(ttInt,ZZy,0,1,PowerSignal);
%%
figure4 = figure('WindowState','maximized');
PowerSignal = 0.999;
ff=get(gca);
     subplot(2,2,1)
     plot(ttInt,ZZx)
     yline(-eta_s,'--k',{'support № 2'},'LineWidth',3,'FontName','Times New Roman','FontSize',14);
     yline(eta_s*0.5,'--k',{'supports № 1,3'},'LineWidth',3,'FontName','Times New Roman','FontSize',14);
     xlabel('\tau','FontName','Times New Roman','FontSize',16)
     ylabel('\xi_{_x}','FontName','Times New Roman','FontSize',16)
     title(['Moves on oX','; ',num2str(nT-tNS),' < nT < ',num2str(nT)],'FontName','Times New Roman','FontSize',16)
     grid on;
     
     subplot(2,2,2)
     h1=stem(fx,Xspectr);
     xline(1/TN,'--r',{'Rotation speed'},'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('f','FontName','Times New Roman','FontSize',16)
     ylabel('Amplitude','FontName','Times New Roman','FontSize',16)
     title(['Spectrum by oX','; ',num2str(PowerSignal*100),'% signal power'],'FontName','Times New Roman','FontSize',16)
     grid on;
     set(get(h1,'BaseLine'),'LineStyle','-');
     set(h1,'MarkerFaceColor','red');
     xlim padded
     ylim padded
     
    
     subplot(2,2,3)
     plot(ttInt,ZZy)
     xlabel('\tau','FontName','Times New Roman','FontSize',16)
     ylabel('\xi_{_y}','FontName','Times New Roman','FontSize',16)
     title(['Moves on oY','; ',num2str(nT-tNS),' < nT < ',num2str(nT)],'FontName','Times New Roman','FontSize',16)
     grid on;
     
     subplot(2,2,4)
     h2=stem(fy,Yspectr);
     xline(1/TN,'--r',{'Rotation speed'},'LineWidth',2,'FontName','Times New Roman','FontSize',14);
     xlabel('f','FontName','Times New Roman','FontSize',16)
     ylabel('Amplitude','FontName','Times New Roman','FontSize',16)
     title(['Spectrum by oY','; ',num2str(PowerSignal*100),'% signal power'],'FontName','Times New Roman','FontSize',16)
     grid on;
     set(get(h2,'BaseLine'),'LineStyle','-');
     set(h2,'MarkerFaceColor','red');
     xlim padded
     ylim padded
     
     %% Запись сил в упругой опоре
 ltt=length(tN); PItog=zeros(2,ltt); flag=1;
 NumbOp=60:120:300;  % Расположение опор (degrees)
 for alfa1=NumbOp
     for k=1:ltt 
        alfa=(alfa1*pi)/180; % переход к радианам
        % Связь перемещений
        uRxs=Z0{IND}(k,mnew+ii+1)*cos(alfa)+Z0{IND}(k,mnew+ii+2)*sin(alfa);
        uRys=-Z0{IND}(k,mnew+ii+1)*sin(alfa)+Z0{IND}(k,mnew+ii+2)*cos(alfa);
        % Связь скоростей
        duRxs=Z0{IND}(4*(mnew+1)+(mnew+1)+ii)*cos(alfa)+Z0{IND}(4*(mnew+1)+(mnew+2)+ii)*sin(alfa);
        duRys=-Z0{IND}(4*(mnew+1)+(mnew+1)+ii)*sin(alfa)+Z0{IND}(4*(mnew+1)+(mnew+2)+ii)*cos(alfa);
        % Связь сил
        Sn(flag,k)=(kappa_s*(uRxs-eta_s)+D_s*duRxs)*(uRxs-eta_s>=0)*((kappa_s*(uRxs-eta_s)+D_s*duRxs)>=0);
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
    subplot(211);hold on;box on;grid on
       for j=1:1:lenOp
           plot(tN,Sn(j,:))
       end
       legend(Legend)
       title('Нормальные реакции в опорах')
       ylabel('S_ { n\alpha_{j}}')
       ff = gca; 
       ff.FontName = 'Times New Roman';
       ff.FontSize = 20;
    
    subplot(212);hold on;box on;grid on
       for j=1:1:lenOp
           plot(tN,St(j,:))
       end
        legend(Legend)
        title('Тангенциальные реакции в опорах')
        xlabel('{\itt} / (2 \pi / N )')
        ylabel('S_ { t\alpha_{j}}')
       ff = gca; 
       ff.FontName = 'Times New Roman';
       ff.FontSize = 20;
     
     