function Right_part = solver(t,Z0,A0,A1,A2,kappa_B, eta_B,amplitude,phase,period)
    
    ksi_0 = amplitude*sin((2*pi/period)*t+phase);
    ksi_0_point = (2*pi/period)*amplitude*cos((2*pi/period)*t+phase);

    MatrKoeff = [zeros(length(A0)) , eye(length(A0)); -A2\A0, -A2\A1]; % матрица коэффициентов
    Q = [zeros((length(A0)-1),1);kappa_B*ksi_0+2*eta_B*ksi_0_point];
    SvobVector = [zeros((length(A0)),1);-A2\Q]; % свободный вектор
    
    Right_part = MatrKoeff*Z0+SvobVector;
end