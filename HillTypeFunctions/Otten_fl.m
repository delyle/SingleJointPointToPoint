function f = Otten_fl(L,b)
% this function computes the normalized isometric contractile force for a given
% length according to a Hill-type model parameterized by Otten 1985
%
% f = exp([ -abs( [L^b(2) - 1]/b(3) ).^b(1) ]), L >= 0
% 
% The parameter b follows Murtola 2022
%
% F is the output force per maximum isometric force
%
% L is the muscle length normalized to optimal length
%
% b is a vector of length 3 of shape parameters
%   b(1) is "roundedness"
%   b(2) is "skewness"
%   b(3) is "width"

t1 = L.^b(2) - 1;
t2 = t1/b(3);
t3 = abs(t2);
t4 = t3.^b(1);
t5 = -t4;
f = exp(t5);