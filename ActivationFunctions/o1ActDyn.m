function da = o1ActDyn(a,u,alp)
% Order-1 activation dynamic model
% Used by Thelan 2003 (https://doi.org/10.1115/1.1531112)
% a is activation
% u is excitation
% alp is a vector of length two
%   alp(1) is the activation time constant (alpha_act in Murtola2023)
%   alp(2) is the deactivation time constant (alpha_deact in Murtola2023)

A_on = alp(1)*(0.5+1.5*a);
A_off = alp(2)./(0.5+1.5*a);

da_on = (u - a)./A_on;
da_off = (u - a)./A_off;

da = if_else(u > a,da_on,da_off);






