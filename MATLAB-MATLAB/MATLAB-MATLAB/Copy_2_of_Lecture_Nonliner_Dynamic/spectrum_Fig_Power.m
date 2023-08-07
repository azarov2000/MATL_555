function [fg,Yg,Trand] = spectrum_Fig_Power(t,x,Ntrand,Fig,eta_Power)
 
 P=polyfit(t,x,Ntrand); Trand=polyval(P,t);
 
 X=x-Trand;
 X=x;
 
 N=length(X); Nf=round(N/2);
 
 dt=mean(diff(t));
 
 % „астота опроса
 fe=1/dt;

 % Ѕыстрое преобразование ‘урье
 fftX=fft(X);     

 % јмплитуды гармоник Y(i) для частот f(i) = i*1/(max(t)-min(t)
 
 A=abs(fftX);
 Y=A(1,1:Nf)/Nf; LY=length(Y);
 % —пектр частот
 f=[0:Nf-1]*fe/N;
 Pow_Y=sum(Y.*Y);
 Pow_eta=0;
 j=0;
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

