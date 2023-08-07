function [DET] = Func(data)
A0 = kron(data.I,(data.E-2*data.zeta_V*data.N*data.R))+data.Mom*kron(data.G03Int,data.R);
DET = det(A0);
end

