function F = Murtola_fp_smooth(L,c,s)
% this function computes the normalized passive tension force
% according to a Hill-type model used by Murtola and Richards 2023
%
% F = c(1)*(exp(c(2)*(L-c(3))-1)    , c(3) < L
%   = 0                             , L <= c(3)
%
%
% F is the output force per maximum isometric force
% L is muscle length divided by optimal length
%
% c is a vector of length 3 of shape parameters
%   c(1) is a scaling parameter
%   c(2) is a rate constant
%   c(3) is the slack length
%
% This function requires CasADi. For an equivalent function without CasADi,
% use Mutola_fp_nc
sigmoid = @(x) 1./(1 + exp(-s*x));

r1 = L-c(3);
r2 = c(2)*r1;
r3 = exp(r2);
r4 = r3-1;
Fp = c(1)*r4;

sigm = sigmoid(L - c(3));
F = Fp.*sigm;


