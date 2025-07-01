function da = o1ActDynSmooth(a,u,alp,s)
% Order-1 activation dynamic model
% Used by Murtola and Richards 2022
% a is activation
% u is excitation
% alp is a vector of length two
%   alp(1) is the activation time constant (alpha_act in Murtola2023)
%   alp(2) is the deactivation time constant (alpha_deact in Murtola2023)
%   s is a smoothing parameter



A_on = alp(1)*(0.5+1.5*a);
A_off = alp(2)./(0.5+1.5*a);

da_on = (u - a)./A_on;
da_off = (u - a)./A_off;

sigmoid = @(x) 1./(1 + exp(-s*x));

da = da_on.*sigmoid(u-a)+da_off.*sigmoid(a-u);

%if_else(u > a,da_on,da_off);




