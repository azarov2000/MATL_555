function disp_param(d)

% Задайте параметры
param1 = 42;
param2 = 'Hello, World!';
param3 = [1, 2, 3];

% Создайте строку с информацией о параметрах
infoString = sprintf('Параметр 1: %d\nПараметр 2: %s\nПараметр 3: %s', param1, param2, mat2str(param3));

% Задайте размер шрифта
fontSize = 14;

% Создайте структуру с настройками окна
h = msgbox(infoString, 'Информация о параметрах', 'help');

% Получите текущий размер окна
currentPosition = get(h, 'Position');

% Увеличьте ширину окна (например, на 100 пикселей)
newWidth = currentPosition(3) + 100;

% Установите новый размер окна
set(h, 'Position', [currentPosition(1), currentPosition(2), newWidth, currentPosition(4)]);

% Получите доступ к объекту текста в окне и установите размер шрифта
txt = findobj(h, 'Type', 'Text');
set(txt, 'FontSize', fontSize);


end

