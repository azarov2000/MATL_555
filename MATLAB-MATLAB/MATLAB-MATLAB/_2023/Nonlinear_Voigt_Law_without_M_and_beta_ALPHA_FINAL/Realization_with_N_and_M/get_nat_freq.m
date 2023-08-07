function [natural_freq] = get_nat_freq(d,Number_freq)

    d.zeta_V = 0;
    d.zeta_e = 0;
    d.N = 0;
    d.N_z = 0;
    d.M_z = 0;

    A0 = kron(d.I,(d.E-2*d.zeta_V*d.N*d.R)) - d.N_z*kron(d.G02Int,d.E) + d.M_z*kron(d.G01Int,d.R)*d.M;

    A1 = 2*d.zeta_V*kron(d.I,d.E)+2*d.zeta_e*kron(d.G00Int,d.E);

    A2 = d.mu_R*kron(d.G00Int,d.E);

    Eig_values = polyeig(A0,A1,A2);
   

    natural_freq_all = sort(imag(Eig_values(find(imag(Eig_values)>=0))));
    
    natural_freq_all_unique = unique(natural_freq_all);

    for i=1:Number_freq
        natural_freq(i) = natural_freq_all_unique(i);
    end 

end

