function [A0,A1,A2,G,D]=get_matrix(N,n,beta,zeta_V,zeta_e)
    
    % Root search
    roots_alpha = solver_for_alpha(1,1,30);

    % Matrix
    E = eye(2); I = eye(n); S=[0 1;-1 0];
    alpha_vect = roots_alpha(1:n);
    for i=1:n
        ALPHA(i,i) = ((alpha_vect(i)))^4;
    end

    % Obtaining orthonormal psi functions
    inc_z = 1000;
    z = linspace(0,1,inc_z);
    dz = mean(diff(z));
    for i=1:n   % index alpha
        K1{i} = @(z) krylovF(z,alpha_vect(i),1);
        K2{i} = @(z) krylovF(z,alpha_vect(i),2);
        K3{i} = @(z) krylovF(z,alpha_vect(i),3);
        K4{i} = @(z) krylovF(z,alpha_vect(i),4);
        Psi_tmp{i} = @(z) (K3{i}(z)-(K3{i}(1)/K4{i}(1))*K4{i}(z));  
        % Normalization
        Int = 0; 
        for j=1:inc_z-1
            dS = 0.5*((Psi_tmp{1}(z(j)))^2+(Psi_tmp{1}(z(j+1)))^2)*dz;
            Int = Int + dS;
        end
        Norm_C(i) = Int;
        Psi{i} = Psi_tmp{i}(z)/sqrt(Norm_C(i));   
    end
    % 
    for i=1:n
        for j = 2:inc_z-1
            DD_Psi{i}(j-1) = (Psi{i}(j-1)+Psi{i}(j+1)-2*Psi{i}(j))/(dz)^2;
        end
        DD_test{i} = @(z) (Norm_C(i))^(-0.5)*(alpha_vect(i))^2*(K1{i}(z)-(K3{i}(1)/K4{i}(1))*K2{i}(z));
    end

    % Filling in the matrix G

    for i=1:n
    for k=1:n
        int_tmp = 0;
        for j=1:inc_z-1
            dS = 0.5*(DD_test{k}(z(j))*Psi{i}(j)+DD_test{k}(z(j+1))*Psi{i}(j+1))*dz;
            int_tmp = int_tmp+dS; 
        end
        G(i,k) = int_tmp;
    end
    end

    % All matrix
    D = kron(2*(zeta_V*ALPHA+zeta_e*I),E);
    A0 = kron(ALPHA,E)+2*zeta_V*N*kron(ALPHA,S);
    A1 = -2*beta*N*kron(G,S)+D; %!!! Минусы перед beta
    A2 = -beta*kron(G,E)+kron(I,E); %!!! Минусы перед beta

end