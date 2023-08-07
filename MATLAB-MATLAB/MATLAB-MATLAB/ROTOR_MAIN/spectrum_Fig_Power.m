function [fg,Yg,Trand] = spectrum_Fig_Power(t,x,Ntrand,eta_Power)
 
P=polyfit(t,x,Ntrand); % ���������� ������������ ��������� �������� ������� (�����.����. ���������.)
Trand=polyval(P,t);    % �������� ���������� �������� ����� ��������
 
X=x-Trand;             % �� ���� �� ������� ���������� ������������ �������
N=length(X);           
Nf=round(N/2);
 
dt=mean(diff(t));      % ��� ������ ��� ����� ��� �����
 
 
fe=1/dt; % ������� ������

fftX=fft(X); % ������� �������������� �����     

A=abs(fftX); 
Y=A(1,1:Nf)/Nf;     % �������� ���������
LY=length(Y);

f=[0:Nf-1]*fe/N;    % ������ ������

Pow_Y=sum(Y.*Y);    % ������ �������� ����������� �������
Pow_eta=0;
j=0;

% ����������� �������� � �����
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

