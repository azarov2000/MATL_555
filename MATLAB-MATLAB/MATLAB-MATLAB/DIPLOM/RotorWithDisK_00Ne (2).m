function RotorWithDisK_00Ne
% Parameters of rotor Dynamics
% Initial values
clear all
close all
global beta zeta_e zeta_i Zeta_e Zeta_te z_C N nameS
% S - indicator of rotor fixation type 
%     == 1 - clamped-hinged supported Rotor
%     == 2 - hinged-hinged supported Rotor
%-----------------------------
   beta=0.005;%025;
 zeta_e=0.025;%25;
 zeta_i=0.05;%05;
 Zeta_e=0.025;%025;
Zeta_te=0.05;%025;

   z_C=input('Disk position (0.3 < z_C < 0.7) = ????  ');
   %-------------------------
   ss='Rotor Support: Clamped-Hinged ''1'' | Hinged-Hinged ''2'' | Clamped-Free ''3''  ';
   nameS=input(ss);
   switch nameS
       case 1
          SupportName='Clamped-Hinged supported Rotor';
       case 2
          SupportName='Hinged-Hinged supported Rotor';
       case 3
          SupportName='Clamped-Free';
   end
   
%-------------------------              
   titlE=['\beta = ',          num2str(beta),...
          ', \zeta_{ i} = ',   num2str(zeta_i),...
          ', \zeta_{ e} = ',   num2str(zeta_e),...
          ', Z_{ e} = ',       num2str(Zeta_e),...
          ', Z_{ \theta e} = ',num2str(Zeta_te),...
          ', \zeta_{ C} = ',   num2str(z_C)];
      
%=============================
Nmin=0;%0;
Nmax=30;%10;
  nN=100;
  Nv=linspace(Nmin,Nmax,nN+1); 
 lNv=length(Nv);
 
%=============================
W=[];   S=[];   V=[];
for j=1:lNv
    N=Nv(j);
    [A2,A1,A0]=RotorMatrixN(z_C,N);
    [v,w,s]=polyeig(A0,A1,A2);
    V=[V;v(:)];  W=[W,w(:)];  S=[S,s(:)];
end

%===========================
% Limits for Graphical presentation:
LimN=[Nmin Nmax];

%---------------------------------
f1=figure('Name','Im(lambda) VS Re(lambda)');
hold on;box on; %axis equal;
for j=1:lNv
    if max(real(W(:,j)))>0
        col='r.';ms=8;
    else
        col='k.';ms=4;
    end
    if j==1
       plot(real(W(:,j)),imag(W(:,j)),'bo','MarkerSize',4,'LineWidth',2)
    else
        plot(real(W(:,j)),imag(W(:,j)),col,'MarkerSize',ms)
    end
end
    xl=get(gca,'XLim');xb=0.1;Xg=XtextPosition(xl,xb); 
    yl=get(gca,'YLim');yb=0.9;Yg=YtextPosition(yl,yb);
    plot(xl,[0 0],'-k');plot([0 0],yl,'-k');
    text(Xg,Yg,SupportName)
    xlabel('RE(\lambda)')
    ylabel('IM(\lambda)','Rotation',0)
    title(titlE)
    grid on

%---------------------------------
f2=figure('Name','Im(lambda) VS N');
hold on;box on
for j=1:lNv
    if max(real(W(:,j)))>0
        col='r.';ms=8;
    else
        col='k.';ms=4;
    end
    plot(Nv(j),imag(W(:,j)),col,'MarkerSize',ms)
end
plot(LimN,[0 0],'-k')
plot([Nmin Nmax],[Nmin  Nmax],'-','LineWidth',1,'Color',[0.07 0.62 1])
plot([Nmin Nmax],[Nmin -Nmax],'-','LineWidth',1,'Color',[0.07 0.62 1])
set(gca,'XLim',[Nmin Nmax]);
    xl=get(gca,'XLim');xb=0.1;Xg=XtextPosition(xl,xb); 
    yl=get(gca,'YLim');yb=0.9;Yg=YtextPosition(yl,yb);
    text(Xg,Yg,SupportName)
    text(Xg,0.3,'IM(\lambda) = N')
    xlabel('N')
    ylabel('IM(\lambda)','Rotation',0)
    title(titlE)
    grid on

%---------------------------------
f3=figure('Name','Re(lambda) VS N');
hold on;box on
for j=1:lNv
    if max(real(W(:,j)))>0
        col='r.';ms=8;
    else
        col='k.';ms=4;
    end
    plot(Nv(j),real(W(:,j)),col,'MarkerSize',ms)
end
plot(LimN,[0 0],'-k')
set(gca,'XLim',[Nmin Nmax]);
    xl=get(gca,'XLim');xb=0.1;Xg=XtextPosition(xl,xb); 
    yl=get(gca,'YLim');yb=0.7;Yg=YtextPosition(yl,yb);
    plot(LimN,[0 0],'-k')
    text(Xg,Yg,SupportName)
    xlabel('N')
    ylabel('RE(\lambda)','Rotation',0)
    title(titlE)
    grid on
    

%***********************************************************************


%=========================================
pause(1)
      N=input('N > N_biff = '); % see on fig 2
      
kappa_s=200;        % Support stiffness
    D_s=0.05;       % Support Damping Parameter
    f_s=0.3; 
  eta_s=0.03;       % Distance to support
     nn=1;          % Exponent in the law of nonlinearity
     eD=0.2;        % Ropport of disk diameter to length: D_disk/l_rotor
      e=0.02;       % Disc-rotor eccentricity

     SS=[0 1;-1 0];  ee=eye(2);  oo=zeros(2,2);
  
title2=[
        'N = ',              num2str(N),...
        ', \zeta_{ C} = ',   num2str(z_C),...
        ', \epsilon_{ D} = ',num2str(eD),...
        ', \eta_{ s} = ',    num2str(eta_s),...
        ', \kappa_{ s} = ',  num2str(kappa_s),...
        ', D_{ s} = ',       num2str(D_s),...
        ', f_{ s} = ',       num2str(f_s),...
        ', \epsilon = ',     num2str(e),...
        ];

[A2,A1,A0]=RotorMatrixN(z_C,N);
B1=A2\A1;B0=A2\A0;
[G00,G0s,Gr0,Grs]=GRINandDERIVATIVES(z_C,z_C,nameS);
%---------------------------------
opt=odeset('AbsTol',1e-10,'RelTol',1e-8);
TN=2*pi/N;      % TN - Time of one revolution (TN=2*pi/N)
nT=1000;        % # of time revolution in implementation
tlim=TN*nT;       
T=[0,tlim];
x0=[1;1;0;0;0;0;0;0]*0.02;      % Initial conditions
[t,Y]=ode45(@(t,Y) rGap(t,Y),T,x0,opt);
tN=t/TN;

%---------------------------------
figure;         % Y VS t/TN ...
    subplot(211);hold on;box on
    plot(tN,Y(:,1),tN,Y(:,2))
    legend('\xi_{ x}','\xi_{ y}')
    title(title2)
    
    subplot(212);hold on;box on
    plot(tN,Y(:,3),tN,Y(:,4))
    legend('\theta_{ x}','\theta_{ y}')
    xlabel('{\itt} / (2 \pi / N )')
%--------------------------------- 
 tNS=50;
 I=find(tN>(tN(end)-tNS));
 figure;        % xi_y VS xu_x       
    hold on;box on
    plot(Y(I,1),Y(I,2))
    axis equal
    at=[num2str(tN(I(end))-tNS),' < {\itt} \rm / (2 \pi / N )< ',num2str(tN(I(end)))];
    xlabel(at)
    title('\xi_{ y} VS \xi_{ x}')
    
%---------------------------------
ltt=length(t);  C=zeros(2,ltt);  F=zeros(2,ltt);
for k=1:ltt 
      Vrel=(N*eD/2-Y(k,6))/abs(N*eD/2-Y(k,6));
    C(:,k)=[1;f_s*Vrel];
     fh(k)=-((Y(k,1)-eta_s)>=0)*(kappa_s*(Y(k,1)-eta_s)^nn+D_s*Y(k,5));
    F(:,k)=fh(k)*C(:,k);
end
figure;
    subplot(211);hold on;box on;grid on
    plot(tN,F(1,:))
    ylabel('F_{ s x}')
    title(title2)
    
    subplot(212);hold on;box on;grid on
    plot(tN,F(2,:))
    xlabel('{\itt} / (2 \pi / N )')
    ylabel('F_{ s y}')
        
%---------------------------------
figure;         % F_x VS t/TN
    subplot(211);hold on;box on;grid on
    plot(tN(I),F(1,I))
    ylabel('F_{ s x}')
    title(title2)
    
    subplot(212);hold on;box on;grid on
    plot(tN(I),F(2,I))
    xlabel('{\itt} / (2 \pi / N )')
    ylabel('F_{ s y}')


%=========================================
function S=rGap(t,Y)
    
        c=[cos(N*t);sin(N*t)];
        P=-e*N*(N*ee-2*Zeta_e*SS)*c;
        Vrel=(N*eD/2-Y(6))/abs(N*eD/2-Y(6));
        C=[1;f_s*Vrel];
        fh=-((Y(1)-eta_s)>=0)*(kappa_s*(Y(1)-eta_s)^nn+D_s*Y(5));
        F=fh*C;
    ff=A2\[G00*(P+F);-Gr0*SS*(P+F)];
    a=Y(1:4);b=Y(5:8);
    S=[b;-B1*b-B0*a+ff];
end
        
%=========================================
function Xg=XtextPosition(xl,xb) %#ok<*DEFNU>
    Xg=xl(1)+xb*(xl(2)-xl(1));
end
function Yg=YtextPosition(yl,yb)
    Yg=yl(1)+yb*(yl(2)-yl(1));
end


end
