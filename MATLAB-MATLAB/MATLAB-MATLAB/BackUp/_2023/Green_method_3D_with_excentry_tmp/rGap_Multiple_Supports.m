function Right_part=rGap_Multiple_Supports(t,Z0,A0,A1,A2,g,h_matr,N,eD,mnew,muR,chi_j,kappa_j,d_j,epsilon,gamma,f_j,NumbOp,disp,disp_vel)

PItog=zeros(2,1); % Подготавливаем вектор сил

%% Учёт силы эксцентриситета
c=[cos(N*t);sin(N*t)];
P = epsilon*N*N*c;            
%% Запись сил от опор


for alfa1=NumbOp
    alfa=(alfa1*pi)/180; % Пересчитываем углы в радианы
    % Связь перемещений
    uRxs=Z0(disp(1))*cos(alfa)+Z0(disp(2))*sin(alfa);
    uRys=-Z0(disp(1))*sin(alfa)+Z0(disp(2))*cos(alfa);
    % Связь скоростей
    duRxs=Z0(disp_vel(1))*cos(alfa)+Z0(disp_vel(2))*sin(alfa);
    duRys=-Z0(disp_vel(1))*sin(alfa)+Z0(disp_vel(2))*cos(alfa);
    % Связь сил
    Sn=(kappa_j*(uRxs-chi_j)+d_j*duRxs)*((uRxs-chi_j)>=0)*((kappa_j*(uRxs-chi_j)+d_j*duRxs)>=0); 
    St=f_j*((N*eD/2+duRys)/abs(N*eD/2+duRys))*Sn;
    MatrRot=[cos(alfa), -sin(alfa); sin(alfa),cos(alfa)];
    PItog=PItog+MatrRot*[-Sn;-St];
end
            
GGG=[-gamma;0]; % Учёт силы тяжести диска

xx=ones(2*(mnew+1),1); % физического смысла сам по себе не несёт, нужен для записи итогового вектора сил

Pext =(P+PItog+GGG); % вектор сил, действующий на узел с диском (2x1)
i_tild = [-muR*gamma;0];
PI=kron(xx,Pext);   % Итоговый вектор ПИ

Gamma = kron(xx,i_tild);

MatrKoeff = [zeros(4*(mnew+1),4*(mnew+1)) , eye(4*(mnew+1),4*(mnew+1)); -A2\A0, -A2\A1]; % матрица коэффициентов
SvobVector = [zeros(4*(mnew+1),1);A2\(g*PI+h_matr*Gamma)]; % свободный вектор

Right_part=MatrKoeff*Z0+SvobVector; % Возврат вектора правой части

 end