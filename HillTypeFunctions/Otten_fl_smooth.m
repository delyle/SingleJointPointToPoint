function f = Otten_fl_smooth(L,b,s)
% this function computes the normalized isometric contractile force for a given
% length according to a Hill-type model parameterized by Otten 1985
%
% f = exp([ -abs_smooth( [L^b(2) - 1]/b(3) ).^b(1) ]), L >= 0
% 
% here we smooth the absolute value function with
%   abs_smooth = x.*tanh(x*s)
% where s is a smoothing parameter (larger = more smooth)
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
t3 = t2.*tanh(t2*s);
t4 = t3.^b(1);
t5 = -t4;
f = exp(t5);