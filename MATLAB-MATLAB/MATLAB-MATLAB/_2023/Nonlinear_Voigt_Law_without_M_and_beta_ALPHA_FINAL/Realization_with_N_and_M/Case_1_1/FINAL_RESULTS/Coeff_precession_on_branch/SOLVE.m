function [RES] = SOLVE(d,T,opt)

number_node = round((d.m-1)/2);
[A0,A1,A2] = Filling_matrix_A(d);


tic;
[t,X] = ode23t(@(t,X) solver_(t,X,A0,A1,A2,d.G02Int,d.G00,d.E,d.zeta_VV,d.alpha,d.M,d.N,d.R,d.m,d.Gamma_0,d.natural_freq,d.F_coeff),T,d.x0,opt);
toc;

for k=1:2:length(A0)
    ksi_X{k} = X(:,k);
end
for k=2:2:length(A0)
    ksi_Y{k} = X(:,k);
end
ksi_X(:,2:2:end) = [];
ksi_Y(:,1:2:end) = [];
ksi_X = [zeros(length(ksi_X{1}),1),ksi_X,zeros(length(ksi_X{1}),1)];
ksi_Y = [zeros(length(ksi_Y{1}),1),ksi_Y,zeros(length(ksi_Y{1}),1)];
z_vect = 0:(1/d.m):1;

for k=1:length(z_vect)
    Z_coord{k} = ones(length(ksi_X{1}),1) * z_vect(k);
end


ind_1 = find(t>=T(end)-d.Number_end);
ind_2 = find(t>=T(end)-(d.Base_number+2*d.Number_end) &...
        t<=T(end)-(d.Base_number+d.Number_end));
RES.MAX_AMPL_1 = max(ksi_X{number_node+1}(ind_1));
RES.MAX_AMPL_2 = max(ksi_X{number_node+1}(ind_2));
[RES.t_extr_All,RES.Extr_X_ALL]=ext(t,ksi_X{number_node+1});
[~,Extr_X_1]=ext(t(ind_1),ksi_X{number_node+1}(ind_1));
[~,Extr_X_2]=ext(t(ind_2),ksi_X{number_node+1}(ind_2));
Extr_X_plus_ind_1 = find(Extr_X_1>=0);
Extr_X_plus_ind_2 = find(Extr_X_2>=0);
RES.EXTR_1 = Extr_X_1(Extr_X_plus_ind_1);
RES.EXTR_2 = Extr_X_2(Extr_X_plus_ind_2);
B = mean(RES.EXTR_1);
A = mean(RES.EXTR_2);
End_cond = X(end,:);

RES.Rel_tol_extr = abs(B-A)/((A+B)/2);
RES.End_cond = End_cond;

[~,CoeffPr] = PrecessionCoeff(t(ind_1), ksi_X{number_node+1}(ind_1),ksi_Y{number_node+1}(ind_1),d.N);

RES.CoffPrecession = mean(CoeffPr);


figure;
box on; grid on; hold on;
plot(t,ksi_X{number_node+1})
ff = gca;
title(['\Omega = ',num2str(d.N),' RelTol = ', num2str(RES.Rel_tol_extr)])
ff.FontSize = 14;

function Right_part = solver_(t,X,A0,A1,A2,G02Int,G00,E,zeta_VV,alpha,M,N,R,m,Gamma_0,nat_freq,F_coeff)

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

    G_g_NL = - alpha*kron(G02Int_g_NL,E);
    F_g_NL = ksi;

    % Force
    EC = zeros(m-1);
    node_N = round((m-1)/2);
    EC(node_N,node_N) = 1;

    F_gamma = -Gamma_0 * cos(F_coeff*nat_freq(1)*t) * [1;0];
    Gamma = kron(G00*EC,E) * kron(ones(m-1,1),F_gamma);


    MatrKoeff = [zeros(length(A0)) , eye(length(A0)); -A2\A0, -A2\A1];
    SvobVector = [zeros((length(A0)),1);A2\(G_f_NL*F_f_NL+G_g_NL*F_g_NL+Gamma)]; 

    Right_part = MatrKoeff*X+SvobVector;
end

end

