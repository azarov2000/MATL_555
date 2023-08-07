close all
clc
clear

m = 6;
h = 1/m;

M = zeros(m-1);

for i=1:length(M)
    for j=1:length(M)
        
        if i==j
            M(i,j) = -2;
        end
        if i == j-1 || j == i-1 
            M(i,j) = 1;
        end 
        
    end 
end

M = (1/h^2)*M;
 
% M0 = [eye(length(a_sub)), zeros(length(a_sub));...
%       zeros(length(a_sub)), zeros(length(a_sub))];
% 
% M1 = [zeros(length(a_sub)), zeros(length(a_sub));...
%       zeros(length(a_sub)), eye(length(a_sub))];
  
% A0 = M0*M;
% A1 = M*M1;
  


  
   