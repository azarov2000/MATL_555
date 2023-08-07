clc
clear
close all
%%

k = 10000;
m_vector = 0:0.0001:10;


i=1;
for m=m_vector
    A = [1, 0,  0,  1;
         0, 1,  0,  1i*m;
         -k,    -k, -k-2*1i*m,  -k*exp(1i*m);
         0, 1,  2,  1i*m*exp(1i*m)];
    x(i) = real(det(A));
    y(i) = imag(det(A));
    i=i+1;
end
%%
f_x = @(xx) interp1(m_vector,x,xx);
f_y = @(xx) interp1(m_vector,y,xx);

fzero(f_x,6.5)
fzero(f_y,6.5)


%%
figure;
hold on; box on; grid on;
plot(m_vector,x,'LineWidth',2)
plot(m_vector,y,'LineWidth',2)
xlabel('M')
ylabel('M,   det(A)')
legend('Re {det(A)}','Im {det(A)}')
title(['\kappa = ', num2str(k)])
ff = gca;
ff.FontSize = 40;

 