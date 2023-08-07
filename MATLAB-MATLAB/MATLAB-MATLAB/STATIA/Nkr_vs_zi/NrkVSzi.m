function [KK]=NrkVSzi(System_name,zi)
global Nmass LenN
WW={};
for j=1:LenN
    N=Nmass(j);U=[];w=[];
     [w] = MatrixOfGreen_Var_N(System_name,N,zi);
    J=find(abs(w)~=inf);U=w(J);
    WW{j}=U(:);
end
lenU=length(U);
%%
flag=0;
KK=0;
for j=1:LenN
   for i=1:lenU
       if real(WW{j}(i))>0
          flag=1;
          KK=j;
          break;
       else 
          flag=0; 
          continue;
       end     
   end
   if flag==1
      break
   end
end



