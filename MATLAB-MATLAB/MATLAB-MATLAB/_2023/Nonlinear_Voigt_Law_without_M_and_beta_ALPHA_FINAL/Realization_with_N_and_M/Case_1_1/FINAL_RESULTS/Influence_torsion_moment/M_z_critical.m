clc 
clear
close all

m_vector = 0:0.0001:3;


i=1;
for m=m_vector
    A = [1, 0,  0,  1;
         0, 1,  0,  1i*m;
         1, 1,  1,  exp(1i*m);
         0, 1,  2,  1i*m*exp(1i*m)];
    x(i) = real(det(A));
    y(i) = imag(det(A));
    i=i+1;
end
%%
f_x = @(xx) interp1(m_vector,x,xx);
f_y = @(xx) interp1(m_vector,y,xx);

fzero(f_x,0.1)
fzero(f_y,0.1)


%%
figure;
hold on; box on; grid on;
plot(m_vector,x,m_vector,y)
xlabel('M')
ylabel('Re [\xi_z],   Im [\xi_z]')
legend('Re [\xi_z]','Im [\xi_z]')

ff = gca;
ff.FontSize = 20;