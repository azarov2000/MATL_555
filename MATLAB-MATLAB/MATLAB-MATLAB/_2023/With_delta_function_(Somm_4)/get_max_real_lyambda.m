function [marker_style,marker_size] = get_max_real_lyambda(data)

a_ksi_0 = kron(data.I,(data.E+2*data.zeta_V*data.N*data.S))+data.kappa_B*kron(data.G00*data.EC,data.E);
a_theta_0 = data.Mom*kron(data.G01*data.EC,data.E)-data.Mom*kron(data.G02Int,data.E);

b_ksi_0 = -data.kappa_B*kron(data.G10*data.EC,data.S);
b_theta_0 = kron(data.I,(data.E+2*data.zeta_V*data.N*data.S))-data.Mom*kron(data.G11*data.EC,data.S)+data.Mom*kron(data.G12Int,data.S);

A0 = [a_ksi_0, a_theta_0;
      b_ksi_0, b_theta_0];

a_ksi_1 = 2*data.zeta_V*kron(data.I,data.E)+2*data.zeta_e*kron(data.G00Int,data.E)+2*data.eta_B*kron(data.G00*data.EC,data.E);
a_theta_1 = 2*data.beta*data.N*kron(data.G00*data.EC,data.E)-2*data.beta*data.N*kron(data.G01Int,data.E);

b_ksi_1 = -2*data.zeta_e*kron(data.G10Int,data.S)-2*data.eta_B*kron(data.G10*data.EC,data.S);
b_theta_1 = 2*data.zeta_V*kron(data.I,data.E)-2*data.beta*data.N*kron(data.G10*data.EC,data.S)+2*data.beta*data.N*kron(data.G11Int,data.S);

A1 = [a_ksi_1, a_theta_1;
      b_ksi_1, b_theta_1];

a_ksi_2 = kron(data.G00Int,data.E)+data.mu_B*kron(data.G00*data.EC,data.E);
a_theta_2 = -data.beta*kron(data.G00*data.EC,data.S)+data.beta*kron(data.G01Int,data.S);

b_ksi_2 = -kron(data.G10Int,data.S)-data.mu_B*kron(data.G10*data.EC,data.S);
b_theta_2 = -data.beta*kron(data.G10*data.EC,data.E)+data.beta*kron(data.G11Int,data.E);

A2 = [a_ksi_2, a_theta_2;
      b_ksi_2, b_theta_2];
  
  
EigenValues = polyeig(A0,A1,A2);


max_real = max(real(EigenValues));

index_number = find(real(EigenValues) == max_real);

marker_style = '.k';
marker_size = 14;
if max_real>0
    marker_style = '.r';
    marker_size = 20;
end 

if max_real>0 && max(imag(EigenValues(index_number))) <= 0.001
    marker_style = '.b';
    marker_size = 20;
end




end

