%________________________
clc; clear
%close all
%___________________ 
Z=zeros(2);E=eye(2);S=[0 1;-1 0];

% stiffness of ...
k = 5;
dk=0;kB=k*[1-dk 0;0 1+dk]*E;

%___Moment limits
M_min=0;    
M_max=10;
dM = 0.0001;
M_vector=M_min:dM:M_max;

%___ initialization
D=[];F=[];

%__ Moments Way:
for M=M_vector
    c=cos(M); s=sin(M);Ex=c*E-s*S;m=M*S;
    DD=[Z, Z, E, E;
        Z, E, Z, -m;
        E, E, Z, -m*Ex;
        m-kB/2, -kB, -kB, -kB*Ex];
    ff=det(DD);
    D=[D,ff];
%    DD
%    pause
end
%% ___ Results prezentation:

figure
    box on;hold on
    plot(M_vector,D,'b-')
    v=gca;
    v.XLim=[M_min;M_max];
    v.YLim=[-1;10];
    plot(v.XLim,[0 0],'k-')
    grid on;
    xlabel('\bf M')
    ylabel('det( {\bfD })')
    title(['\kappa = ',num2str(k) ...
           ':  \kappa_{ x} = ',num2str(k-dk),...
           ',  \kappa_{ y} = ',num2str(k+dk),...
           ',  \Delta_{ \kappa} = ',num2str(dk),...
           ',  \delta_{ M} = ',num2str(dM)])
max(D),min(D)