function Right_part = solver(t,Z0,dat)
    
    % Дисбаланс
    epsilon = @(x) dat.Ampl_eps*sin(pi*x);
    pfi = @(x) dat.Ampl_phase*sin(2*pi*x);
    h = 1/dat.m;
    absciss = h:h:1;
    f = [];
    for i=1:length(absciss)
        f = [f;epsilon(absciss(i))*dat.N*dat.N*cos(dat.N*t+pfi(absciss(i)))];
        f = [f;epsilon(absciss(i))*dat.N*dat.N*sin(dat.N*t+pfi(absciss(i)))];
    end
    
    % Вектор центробежных сил
    Pg = dat.G*kron(ones(2,1),f);
    
    % Вектор кинематического возбуждения
    ksi_0_x = dat.amplitude*sin((2*pi/dat.period)*t+dat.phase);
    ksi_0_x_point = (2*pi/dat.period)*dat.amplitude*cos((2*pi/dat.period)*t+dat.phase);
    ksi_0_y = dat.amplitude*cos((2*pi/dat.period)*t+dat.phase);
    ksi_0_y_point = -(2*pi/dat.period)*dat.amplitude*sin((2*pi/dat.period)*t+dat.phase);
    pQ = [zeros((length(dat.A0)-2),1);dat.kappa_B*ksi_0_x+2*dat.eta_B*ksi_0_x_point;dat.kappa_B*ksi_0_y+2*dat.eta_B*ksi_0_y_point];
    
    % Итоговый вектор сил (вектор в правой части)
    Q = Pg+pQ;
    
    MatrKoeff = [zeros(length(dat.A0)) , eye(length(dat.A0)); -dat.A2\dat.A0, -dat.A2\dat.A1]; % матрица коэффициентов
    SvobVector = [zeros((length(dat.A0)),1);-dat.A2\Q]; % свободный вектор
    
    Right_part = MatrKoeff*Z0+SvobVector;
end