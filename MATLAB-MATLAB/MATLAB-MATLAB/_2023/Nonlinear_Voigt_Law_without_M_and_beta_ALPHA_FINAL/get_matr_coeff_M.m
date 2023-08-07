function M = get_matr_coeff_M(m)
h = 1/m;
M = zeros(m-1);

for i=1:length(M)
    for j=1:length(M)
        
        if i==j
            M(i,j) = -2;
        end
        if i == j-1 || j == i-1 
            M(i,j) = 1;
        end 
        
    end 
end
M = (1/h)^2*M;
M = kron(M,eye(2));
end

