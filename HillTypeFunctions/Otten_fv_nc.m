function F = Otten_fv_nc(V,d)
% this function computes the normalized contractile force at optimal length
% for a given contractile velocity according to a Hill-type model
% parameterized by Otten 1985
%
% NOTE that contractile velocity is shortening velocity (positive is
% shortening)
%
% f = d(2) - (d(2)-1)*(1+V)./(1-d(3)*V) ,      V < 0
%   = 1 - V./(1+d(1)*V)                 , 0 <= V < 1
%   = 0                                 , 1 <  V
%
% The parameter d follows Murtola 2022
%
% F is the output force per maximum isometric force
%
% V is -(dl/dt)/l_0/V_max, where
%   l is muscle length
%   l_0 is optimal length
%   V_max is maximum contractile velocity
%
% d is a vector of length 3 of shape parameters
%   d(1) is rate of decay in contraction
%   d(2) is maximum extension force
%   d(3) is saturation rate in extension
%
% This function will not behave correctly with CasADi.
% For a function using CasADi, refer to Otten_fv

F = NaN(size(V));

for i = 1:numel(V)
    Vi = V(i);
    if Vi >= 1
        F(i) = 0;
    elseif Vi < 0
        s1 = 1 + Vi;
        s2 = 1 - d(3)*Vi;
        s3 = d(2)-1;
        s4 = s3*s1./s2;
        F(i) = d(2) - s4;
    else
        t1 = 1 + d(1)*Vi;
        t2 = 1 - Vi;
        F(i) = t2./t1;
    end
end
