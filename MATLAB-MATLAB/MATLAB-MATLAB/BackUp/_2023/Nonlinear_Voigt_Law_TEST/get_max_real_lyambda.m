function [marker_style,marker_size,res] = get_max_real_lyambda(data)

A0 = kron(data.I,(data.E-2*data.zeta_V*data.N*data.R))+data.Mom*kron(data.G03Int,data.R);

A1 = 2*data.zeta_V*kron(data.I,data.E)+2*data.beta_R*data.N*kron(data.G02Int,data.R)+2*data.zeta_e*kron(data.G00Int,data.E);

A2 = data.mu_R*kron(data.G00Int,data.E)-data.beta_R*kron(data.G02Int,data.E);

  
  
EigenValues = polyeig(A0,A1,A2);


max_real = max(real(EigenValues));

index_number = find(real(EigenValues) == max_real);

res = EigenValues(index_number);

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

