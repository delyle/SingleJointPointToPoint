function F = Otten_fv_smooth(V,d,s)
% this function computes the normalized contractile force at optimal length
% for a given contractile velocity according to a Hill-type model
% parameterized by Otten 1985
%
% In this version, discontinuities are smoothed using sigmoid functions,
% and singularities are canceled by smooth ramp functions.
%
% The parameter s > 0 is applied to all sigmoid functions (min. 200 recommended)
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


sigmoid = @(x) 1./(1 + exp(-s*x));

sigmV = sigmoid(-V);
sigpV = sigmoid(V);

smoothRampNeg = +1*d(3)*V.*sigpV;
s1 = 1 + V;
s2 = 1 - d(3)*V + smoothRampNeg;
s3 = d(2)-1;
s4 = s3*s1./s2;
Fv_neg = d(2) - s4;

smoothRampPos = -1*d(1)*V.*sigmV;
t1 = 1 + d(1)*V + smoothRampPos;
t2 = 1 - V;
Fv_pos = t2./t1;

F = (Fv_pos.*sigpV + Fv_neg.*sigmV).*sigmoid(1-V);
