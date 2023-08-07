     function [tPr,CoeffPr] = PrecessionCoeff(t_init, x_init,y_init,N)
           length_t=length(t_init);
           tInit=linspace(t_init(1),t_init(end),length_t); % интерполируем с равным шагом
           xInit=interp1(t_init,x_init,tInit);
           yInit=interp1(t_init,y_init,tInit);

            h=(tInit(end)-tInit(1))/(length_t-1); % шаг, с которым интерполировали
            tPr=tInit(2:(end-1));
              X=xInit(2:(end-1));
              Y=yInit(2:(end-1));
              D=X.^2+Y.^2;
                DX=xInit(3:end)-xInit(1:(end-2));
                DY=yInit(3:end)-yInit(1:(end-2));           
                fr=(DY.*X-DX.*Y)./D;
                CoeffPr=fr/(2*h*N);

     end