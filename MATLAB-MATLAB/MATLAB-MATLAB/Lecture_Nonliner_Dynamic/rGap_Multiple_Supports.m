function Right_part=rGap_Multiple_Supports(t,Z0,A0,A1,A2,g,N,eD,mnew,eta_s,kappa_s,nn,D_s,e,gd,f_s,NumbOp,System_name,disp,angle,disp_vel,angle_vel)

%% Учёт силы эксцентриситета
   c=[cos(N*t);sin(N*t)];
   P=e*N*N*c;            
%% Запись сил от опор

PItog=zeros(2,1); % Задание размерности для итогового вектора сил

    switch System_name
    
        case 4 % Заделка – заделка
            
            for alfa1=NumbOp
                alfa=(alfa1*pi)/180; % Пересчитываем углы в радианы
                % Связь перемещений
                uRxs=Z0(disp(1))*cos(alfa)+Z0(disp(2))*sin(alfa);
                uRys=-Z0(disp(1))*sin(alfa)+Z0(disp(2))*cos(alfa);
                % Связь скоростей
                duRxs=Z0(disp_vel(1))*cos(alfa)+Z0(disp_vel(2))*sin(alfa);
                duRys=-Z0(disp_vel(1))*sin(alfa)+Z0(disp_vel(2))*cos(alfa);
                % Связь сил
                Sn=(kappa_s*(uRxs-eta_s)+D_s*duRxs)*((uRxs-eta_s)>=0)*((kappa_s*(uRxs-eta_s)+D_s*duRxs)>=0); 
                St=f_s*((N*eD/2+duRys)/abs(N*eD/2+duRys))*Sn;
                MatrRot=[cos(alfa), -sin(alfa); sin(alfa),cos(alfa)];
                PItog=PItog+MatrRot*[-Sn;-St];
            end

        case 1 % Заделка – шарнир
            
            for alfa1=NumbOp
                alfa=(alfa1*pi)/180;
                % Связь перемещений
                uRxs=Z0(mnew+ii)*cos(alfa)+Z0(mnew+ii+1)*sin(alfa);
                uRys=-Z0(mnew+ii)*sin(alfa)+Z0(mnew+ii+1)*cos(alfa);
                % Связь скоростей
                duRxs=Z0(4*(mnew+1)+(mnew)+ii)*cos(alfa)+Z0(4*(mnew+1)+(mnew+1)+ii)*sin(alfa);
                duRys=-Z0(4*(mnew+1)+(mnew)+ii)*sin(alfa)+Z0(4*(mnew+1)+(mnew+1)+ii)*cos(alfa);
                % Связь сил
                Sn=(kappa_s*(uRxs-eta_s)+D_s*duRxs)*((uRxs-eta_s)>=0)*((kappa_s*(uRxs-eta_s)+D_s*duRxs)>=0); 
                St=f_s*((N*eD/2+duRys)/abs(N*eD/2+duRys))*Sn;
                MatrRot=[cos(alfa), -sin(alfa); sin(alfa),cos(alfa)];
                PItog=PItog+MatrRot*[-Sn;-St];
            end
            
            case 3 % Консольное закрепление
            
            for alfa1=NumbOp
                alfa=(alfa1*pi)/180;
                % Связь перемещений
                uRxs=Z0(mnew+ii)*cos(alfa)+Z0(mnew+ii+1)*sin(alfa);
                uRys=-Z0(mnew+ii)*sin(alfa)+Z0(mnew+ii+1)*cos(alfa);
                % Связь скоростей
                duRxs=Z0(4*(mnew+1)+(mnew)+ii)*cos(alfa)+Z0(4*(mnew+1)+(mnew+1)+ii)*sin(alfa);
                duRys=-Z0(4*(mnew+1)+(mnew)+ii)*sin(alfa)+Z0(4*(mnew+1)+(mnew+1)+ii)*cos(alfa);
                % Связь сил4*(mnew+1)+(mnew)+ii
                Sn=(kappa_s*(uRxs-eta_s)+D_s*duRxs)*((uRxs-eta_s)>=0)*((kappa_s*(uRxs-eta_s)+D_s*duRxs)>=0); 
                St=f_s*((N*eD/2+duRys)/abs(N*eD/2+duRys))*Sn;
                MatrRot=[cos(alfa), -sin(alfa); sin(alfa),cos(alfa)];
                PItog=PItog+MatrRot*[-Sn;-St];
            end
    end
            
  
GGG=[-gd;0]; % Учёт силы тяжести диска

xx=ones(2*(mnew+1),1); % физического смысла сам по себе не несёт, нужен для записи итогового вектора сил

Pext=(P+PItog+GGG); % вектор сил, действующий на узел с диском 2x1
PI=kron(xx,Pext);   % Итоговый вектор P с тильдой

MatrKoeff=[zeros(4*(mnew+1),4*(mnew+1)) , eye(4*(mnew+1),4*(mnew+1)); -A2\A0, -A2\A1]; % матрица коэффициентов
SvobVector=[zeros(4*(mnew+1),1);(A2\g)*PI]; % свободный вектор

Right_part=MatrKoeff*Z0+SvobVector; % Возврат вектора правой части

 end