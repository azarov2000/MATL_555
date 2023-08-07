function data = update_data_structure(I, E, S, zeta_V, Mom, kappa_B, zeta_e, G00Int, ...
    beta, G02Int, eta_B, mu_B, EC, G00, G12Int, G10, G01Int, G11Int, G10Int, G01, G11)
    % Создание структуры
    data = struct();

    % Запись данных в структуру
    data.I = I;
    data.E = E;
    data.S = S;
    data.zeta_V = zeta_V;
    data.Mom = Mom;
    data.kappa_B = kappa_B;
    data.zeta_e = zeta_e;
    data.G00Int = G00Int;
    data.beta = beta;
    data.G02Int = G02Int;
    data.eta_B = eta_B;
    data.mu_B = mu_B;
    data.EC = EC;
    data.G00 = G00;
    data.G12Int = G12Int;
    data.G10 = G10;
    data.G01Int = G01Int;
    data.G11Int = G11Int;
    data.G10Int = G10Int;
    data.G01 = G01;
    data.G11 = G11;
end



