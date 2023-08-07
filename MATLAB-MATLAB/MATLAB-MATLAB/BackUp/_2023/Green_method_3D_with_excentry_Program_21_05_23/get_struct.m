function [dat,sol] = get_struct(m,N,Mom,d,l,ro,Elastic,zeta_e,zeta_V,m_B,k_B,eta_B,...
                            Ampl_eps, Ampl_phase,...
                            amplitude, period, phase,nT,initial_vect)
                       
% Структура dat.
J = pi*d^4/64;          %[m^4]
C = Elastic*J;          %
eps_d = d/l;            %[-]
beta = (eps_d^2)/16;    %[-]
epsA = Ampl_eps/l;

if ~isempty(k_B{1}) && isempty(k_B{2})
    kappa_B = k_B{1}*l^3/C;
elseif isempty(k_B{1}) && ~isempty(k_B{2})
    kappa_B = k_B{2}*48;
else 
    error('Некорректно задана жесткость')
end 

dat = struct(...
    'm', m, ...
    'N', N, ...
    'Mom', Mom, ...
    'zeta_e', zeta_e, ...
    'zeta_V', zeta_V, ...
    'beta', beta, ...
    'mu_B', m_B/(ro*(pi*d^2/4)*l*l), ...
    'kappa_B', kappa_B, ...
    'eta_B', eta_B, ...
    'Ampl_eps', epsA, ...
    'Ampl_phase',Ampl_phase ...
);
dat.amplitude = amplitude;
dat.period = period;
dat.phase = phase;

% Структура sol.
TN=2*pi/N;
tlim=TN*nT;

sol.TN = TN;
sol.nT = nT;
sol.T = [0,tlim];
if isempty(initial_vect)
    sol.x0 = zeros(((2*m+2)*2),1);
else
    sol.x0 = x0;
end