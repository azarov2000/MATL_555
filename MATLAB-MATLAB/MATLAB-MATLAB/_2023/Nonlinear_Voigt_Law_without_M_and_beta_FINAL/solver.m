function Right_part = solver(t,X,A0,A1,A2,G02Int,E,zeta_VV,M,N,R,m)

%% Filling in vectors
ksi = M*X(1:1:2*(m-1));
ksi_dot = M*X(2*(m-1)+1:1:4*(m-1));
t
for i=1:m-1
    ksi_local = ksi(2*i-1:2*i);
    ksi_dot_local = ksi_dot(2*i-1:2*i);
    f_NL(i) = (ksi_dot_local)' * ksi_dot_local - 2*N*(ksi_dot_local)' * R * ksi_local + N^2 * (ksi_local)' * ksi_local;
    
    
    G02Int(:,i) = G02Int(:,i) * f_NL(i);
end
G = 2*zeta_VV*kron(G02Int,E);
F = (-ksi_dot + kron(eye(m-1),R) * N * ksi);



MatrKoeff = [zeros(length(A0)) , eye(length(A0)); -A2\A0, -A2\A1]; % matrix of coefficients
SvobVector = [zeros((length(A0)),1);A2\(G*F)];                     % free vector

Right_part = MatrKoeff*X+SvobVector;
end