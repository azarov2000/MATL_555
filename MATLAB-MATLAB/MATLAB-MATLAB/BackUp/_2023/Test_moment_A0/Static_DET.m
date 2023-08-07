clc
clear
close all

vector_M = linspace(-20,20,1000);

for i=1:length(vector_M)
    DET(i) = get_det(vector_M(i));
end 

figure;
box on; grid on; hold on;
plot(vector_M,real(DET),'LineWidth',1.5)
plot(vector_M,imag(DET),'LineWidth',1.5)
yline(0,'LineWidth',1)
ff = gca;
ff.FontSize = 25;
xlabel('M')
ylabel('det')
legend('Re(det)','Im(det)')


function DET = get_det(M)
    A = [1, 0,  0,  1;
         0, 1,  0,  1i*M;
         1, 1,  1,  exp(1i*M);
         0, 1,  2,  1i*M*exp(1i*M)];
     DET = det(A);
end