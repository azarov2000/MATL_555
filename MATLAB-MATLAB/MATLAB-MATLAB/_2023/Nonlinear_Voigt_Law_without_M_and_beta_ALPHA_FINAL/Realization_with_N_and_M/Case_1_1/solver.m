function Right_part = solver(t,X,A0,A1,A2,G02Int,G00,E,zeta_VV,alpha,M,N,R,m,Gamma_0,nat_freq)

%% Filling in vectors
ksi = M*X(1:1:2*(m-1));
ksi_dot = M*X(2*(m-1)+1:1:4*(m-1));
for i=1:m-1
    ksi_local = ksi((2*i-1):2*i);
    ksi_dot_local = ksi_dot((2*i-1):2*i);
    f_NL(i) = (ksi_dot_local)' * ksi_dot_local - 2*N*(ksi_dot_local)' * R * ksi_local + N^2 * (ksi_local)' * ksi_local;
    g_NL(i) = (ksi_local)' * ksi_local;
    
    G02Int_f_NL(:,i) = G02Int(:,i) * f_NL(i);
    G02Int_g_NL(:,i) = G02Int(:,i) * g_NL(i);
end
G_f_NL = 2*zeta_VV*kron(G02Int_f_NL,E);
F_f_NL = (-ksi_dot + kron(eye(m-1),R) * N * ksi);

% G_g_NL = - alpha*kron(G02Int_g_NL,R);
G_g_NL = - alpha*kron(G02Int_g_NL,E);
F_g_NL = ksi;

% Force
EC = zeros(m-1);
node_N = round((m-1)/2);
EC(node_N,node_N) = 1;

F_gamma = -Gamma_0 * cos(0.8*nat_freq(1)*t) * [1;0];
Gamma = kron(G00*EC,E) * kron(ones(m-1,1),F_gamma);

% F_gamma = -Gamma_0 * cos(20*t) * [0;0;0;0;1;0;0;0;0;0];
% Gamma = kron(G00*EC,E) * F_gamma;
% Gamma = kron(G00*EC,E) * kron(ones(length(A0),1),-Gamma_0 * cos(20*t));


MatrKoeff = [zeros(length(A0)) , eye(length(A0)); -A2\A0, -A2\A1]; % matrix of coefficients
% SvobVector = [zeros((length(A0)),1);A2\(G_f_NL*F_f_NL+G_g_NL*F_g_NL)];                     % free vector
SvobVector = [zeros((length(A0)),1);A2\(G_f_NL*F_f_NL+G_g_NL*F_g_NL+Gamma)];                 % free vector

Right_part = MatrKoeff*X+SvobVector;
end