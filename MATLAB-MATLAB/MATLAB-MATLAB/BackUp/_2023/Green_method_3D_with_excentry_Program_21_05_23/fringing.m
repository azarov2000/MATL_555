function [dat] = fringing(dat)

L = dat.m ;
dat.F0 = dat.F0(2:L);
dat.F2 = dat.F2(2:L);
dat.Fg = dat.Fg(2:L);
dat.T = dat.T(2:L);
dat.G00Int = dat.G00Int(2:L,2:L);
dat.G02Int = dat.G02Int(2:L,2:L);
dat.G03Int = dat.G03Int(2:L,2:L);
dat.G30Int = dat.G30Int(2:L,2:L);
dat.I = dat.I(2:L,2:L);
dat.Z = dat.Z(2:L);

end

