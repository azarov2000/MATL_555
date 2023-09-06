function [natural_freq] = get_nat_freq(d,Number_freq)

d.zeta_V = 0;
d.zeta_e = 0;
d.N = 0;
d.N_z = 0;
d.M_z = 0;

[A0,A1,A2] = Filling_matrix_A(d);

Eig_values = polyeig(A0,A1,A2);
   

natural_freq_all = sort(imag(Eig_values(find(imag(Eig_values)>=0))));

natural_freq_all_unique = natural_freq_all(1);
for i=2:length(natural_freq_all)
   if abs(natural_freq_all(i)-natural_freq_all(i-1))>10e-5
       natural_freq_all_unique = [natural_freq_all_unique,natural_freq_all(i)];
   end
end

for i=1:Number_freq
    natural_freq(i) = natural_freq_all_unique(i);
end 

end

