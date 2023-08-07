function [data] = Filling_structure(I,E,R,zeta_V,zeta_VV,zeta_e,G00Int,G02Int,mu_R,alpha,M,m)

    data.I = I;
    data.E = E; 
    data.R = R;
    data.zeta_V = zeta_V;
    data.zeta_e = zeta_e;
    data.G00Int = G00Int;
    data.G02Int = G02Int;
    data.mu_R = mu_R;
    data.alpha = alpha;
    data.M = M;
    data.m = m;
    data.zeta_VV = zeta_VV;
    
end

