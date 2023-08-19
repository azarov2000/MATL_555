function status = odeoutput(t,y, flag, tstart, tspan)
    persistent lastUpdate
    persistent h
    if nargin < 3
        flag = '';
    end

    switch flag
        case 'init'
            lastUpdate = tic;
            h = waitbar(0, 'Выполняется расчет...');
        case 'done'
            close(h);
            return;
        otherwise
            % Обновление waitbar каждые 0.5 секунды
            if toc(lastUpdate) >= 0.5
                waitbar((t - tstart) / tspan, h);
                lastUpdate = tic;
            end
    end

    status = 0;
end
