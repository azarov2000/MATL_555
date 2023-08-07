function K=krylovF(z,a,N)

%Вычилсяется фукциz Крылова #N, имеющей аргумент a*z,
%   где 
%   а - частотный параметр
%   z - координата, изменяющаяся в пределах
%       0<=z<=1

if N>4 | N<1
    error('Wrong number of function')
end
if z<0 | z>1
    error('Wrong for coordinate z')
end

x=a*z;
switch N
    case 1
        K=cosh(x)+cos(x);
    case 2
        K=sinh(x)+sin(x);
    case 3
        K=cosh(x)-cos(x);
    case 4
        K=sinh(x)-sin(x);
end
K=K/2;