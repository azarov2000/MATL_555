function [x0] = Initial_condition(d,Initial_value)

x0 = zeros(2*(2*d.m-2),1);
number_node = round((d.m-1)/2);
x0(2*number_node-1) = Initial_value;
x0(2*number_node) = Initial_value;

end

