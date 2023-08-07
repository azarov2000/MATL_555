function [fg,Yg,Trand] = spectrum_Fig_Power(t,x,Ntrand,eta_Power)
 
P=polyfit(t,x,Ntrand); % Возвращает коэффициенты полнинома заданной степени (метод.наим. квадратов.)
Trand=polyval(P,t);    % Получаем записываем ординаты этого полинома
 
X=x-Trand;             % по сути мы убираем постоянную составляющую сигнала
N=length(X);           
Nf=round(N/2);
 
dt=mean(diff(t));      % вот именно для этого нам нужна
 
 
fe=1/dt; % Частота опроса

fftX=fft(X); % Быстрое преобразование Фурье     

A=abs(fftX); 
Y=A(1,1:Nf)/Nf;     % Получаем амплитуды
LY=length(Y);

f=[0:Nf-1]*fe/N;    % Спектр частот

Pow_Y=sum(Y.*Y);    % Полная мощность полученного сигнала
Pow_eta=0;
j=0;

% Накапливаем мощность в цикле
 while Pow_eta<=eta_Power*Pow_Y
     j=j+1;
     if j==LY
         break
     else
         Pow_eta=Pow_eta+Y(j)*Y(j);
     end
 end
 fg=f(1:j);Yg=Y(1:j);
 
end

