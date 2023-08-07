function [T,X]=ext(t,x)

if length(t)<=2
   T=t;X=x;
else
   s=sign(diff(x));
   ds=s(1:end-1);Ds=s(2:end);
   I=find(ds~=Ds);
   if isempty(I)
       T=t([1,end]);X=x([1,end]);
   else
       T=t(I+1);X=x(I+1);
   end
end
return

