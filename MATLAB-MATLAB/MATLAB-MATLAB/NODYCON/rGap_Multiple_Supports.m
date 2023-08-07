function SSS=rGap_Multiple_Supports(t,Z0,DD,A0,A1,A2,N,eD,mnew,eta_s,kappa_s,nn,D_s,e,gd,f_s,NumbOp,ii,System_name)

%% Учёт силы эксцентриситета
   c=[cos(N*t);sin(N*t)];
   P=e*N*N*c;            
%% Запись сил от опор

PItog=zeros(2,1); % Задание размерности для итогового вектора сил

    switch System_name
    
        case 4 % Заделка – заделка
            
            for alfa1=NumbOp
                alfa=(alfa1*pi)/180;
                % Связь перемещений
                uRxs=Z0(mnew+ii+1)*cos(alfa)+Z0(mnew+ii+2)*sin(alfa);
                uRys=-Z0(mnew+ii+1)*sin(alfa)+Z0(mnew+ii+2)*cos(alfa);
                % Связь скоростей
                duRxs=Z0(4*(mnew+1)+(mnew+1)+ii)*cos(alfa)+Z0(4*(mnew+1)+(mnew+2)+ii)*sin(alfa);
                duRys=-Z0(4*(mnew+1)+(mnew+1)+ii)*sin(alfa)+Z0(4*(mnew+1)+(mnew+2)+ii)*cos(alfa);
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
%                 if alfa1==NumbOp(1)
%                     eta_s = eta_s1 - (gd/3);
%                 end
%                 if alfa1==NumbOp(2)
%                     eta_s = eta_s1 + (gd/3);
%                 end
                
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
            
  
GGG=[-gd;0]; % Учёт силы тяжести

xx=ones(2*(mnew+1),1); % физического смысла сам по себе не несёт, нужен для записи итогового вектора сил

Pext=(P+PItog+GGG); % вектор сил, действующий на узел с диском
PP=kron(xx,Pext);   % Итоговый вектор P с тильдой

B1=A2\A1; B2=A2\A0; TT=(A2\DD)*PP;
MatrKoeff=[zeros(4*(mnew+1),4*(mnew+1)) , eye(4*(mnew+1),4*(mnew+1)); -B2 , -B1]; % матрица коэффициентов
SvobVector=[zeros(4*(mnew+1),1);TT]; % свободный вектор

SSS=MatrKoeff*Z0+SvobVector; % Возврат вектора правой части

 end