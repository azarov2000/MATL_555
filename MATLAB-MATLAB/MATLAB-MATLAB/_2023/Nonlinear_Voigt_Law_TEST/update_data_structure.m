function data = update_data_structure(I, E, R, zeta_V, Mom, zeta_e, G00Int, beta_R, ...
    G02Int, G03Int, mu_R)
    % Создание структуры
    data = struct();

    % Запись данных в структуру
    data.I = I;
    data.E = E;
    data.R = R;
    data.zeta_V = zeta_V;
    data.Mom = Mom;
    data.zeta_e = zeta_e;
    data.G00Int = G00Int;
    data.beta_R = beta_R;
    data.G02Int = G02Int;
    data.G03Int = G03Int;
    data.mu_R = mu_R;
end



