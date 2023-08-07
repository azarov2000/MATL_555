function Right_part = solver(t,Z0,A0,A1,A2,G,h,N)
    
    % Дисбаланс
    epsilon = @(x) (1e-3/0.7)*sin(pi*x);
    pfi = @(x) 1;%*sin(2*pi*x);
    absciss = h:h:1;
    q = [];
    for i=1:length(absciss)-1
        q = [q;epsilon(absciss(i))*N*N*cos(N*t+pfi(absciss(i)))];
        q = [q;epsilon(absciss(i))*N*N*sin(N*t+pfi(absciss(i)))];
    end
    Q = G*q;
    
    MatrKoeff = [zeros(length(A0)) , eye(length(A0)); -A2\A0, -A2\A1]; % матрица коэффициентов
    SvobVector = [zeros((length(A0)),1);-A2\Q]; % свободный вектор
    
    Right_part = MatrKoeff*Z0+SvobVector;
end