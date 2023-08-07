function Right_part = solver(t,Z0,A0,A1,A2,g_psi,g_theta,g_X,h,N,kappa_B,eta_B,amplitude,phase,period)
    
    % Дисбаланс
    epsilon = @(x) (1e-3/0.7)*sin(pi*x);
    pfi = @(x) 1*sin(2*pi*x);
    absciss = h:h:(1-h);
    f = [];
    for i=1:length(absciss)
        f = [f;epsilon(absciss(i))*N*N*cos(N*t+pfi(absciss(i)))];
        f = [f;epsilon(absciss(i))*N*N*sin(N*t+pfi(absciss(i)))];
    end
%     Pg = g*f;
    
    % Кинематическое возбуждение
    ksi_0_x = amplitude*sin((2*pi/period)*t+phase);
    ksi_0_x_point = (2*pi/period)*amplitude*cos((2*pi/period)*t+phase);
    ksi_0_y = amplitude*cos((2*pi/period)*t+phase);
    ksi_0_y_point = -(2*pi/period)*amplitude*sin((2*pi/period)*t+phase);
    p_ex = [kappa_B*ksi_0_x+2*eta_B*ksi_0_x_point;kappa_B*ksi_0_y+2*eta_B*ksi_0_y_point];
    % Итоговый вектор 
    pf = g_X*f + p_ex;
    pG = [g_psi, zeros(length(g_psi)),zeros(length(g_psi),2);
         zeros(length(g_psi)), g_theta, zeros(length(g_psi),2);
         zeros(2,length(g_psi)),zeros(2,length(g_psi)),eye(2)];
     
    Q = pG*[kron(ones(2,1),f);pf];
    
    MatrKoeff = [zeros(length(A0)) , eye(length(A0)); -A2\A0, -A2\A1]; % матрица коэффициентов
    SvobVector = [zeros((length(A0)),1);-A2\Q]; % свободный вектор
    
    Right_part = MatrKoeff*Z0+SvobVector;
end