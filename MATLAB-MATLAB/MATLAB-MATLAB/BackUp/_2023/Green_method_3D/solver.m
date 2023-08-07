function Right_part = solver(t,Z0,A0,A1,A2,kappa_B, eta_B,amplitude,phase,period)
    
    ksi_0_x = amplitude*sin((2*pi/period)*t+phase);
    ksi_0_x_point = (2*pi/period)*amplitude*cos((2*pi/period)*t+phase);
 
    ksi_0_y = amplitude*cos((2*pi/period)*t+phase);
    ksi_0_y_point = -(2*pi/period)*amplitude*sin((2*pi/period)*t+phase);
    
    MatrKoeff = [zeros(length(A0)) , eye(length(A0)); -A2\A0, -A2\A1]; % матрица коэффициентов
    Q = [zeros((length(A0)-2),1);kappa_B*ksi_0_x+2*eta_B*ksi_0_x_point;kappa_B*ksi_0_y+2*eta_B*ksi_0_y_point];
    SvobVector = [zeros((length(A0)),1);-A2\Q]; % свободный вектор
    
    Right_part = MatrKoeff*Z0+SvobVector;
end