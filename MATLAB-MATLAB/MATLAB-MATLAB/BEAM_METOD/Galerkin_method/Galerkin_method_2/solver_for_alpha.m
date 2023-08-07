function rootS = solver_for_alpha(type_BC,start,finish)
    % Search for the alpha parameter by the type of boundary conditions
    K1 = @(a) krylovF(1,a,2);
    K2 = @(a) krylovF(1,a,2);
    K3 = @(a) krylovF(1,a,3);
    K4 = @(a) krylovF(1,a,4);
    
    switch type_BC
        case 1      % Заделка-заделка
            A =  @(a)(det([K3(a),K4(a);K2(a),K3(a)]));
        case 2      % Шарнир-шарнир
            A = @(a)(det([K2(a),K4(a);K4(a),K2(a)]));
        otherwise
            error ('System Pinning type error') 
    end 
    n_points = 10000;                       % кол-во точек для начальных приближений
    inc = linspace(start,finish,n_points);
    for i=1:n_points
        sol(i) = fzero(A,inc(i));
    end
    rootS = uniquetol(sol,1e-4);    
    
end