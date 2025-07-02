function f_absSmooth = absSmoothFun(type)
% returns a function handle for a smoothed absolute value function
% depending on type (case insensitive)
%
% Each returned function can be implemented as f_absSmooth(x,s), 
% where s is a smoothing parameter s.t. f_absSmooth(x,s) -> |x| as s -> Inf
% 
% Types are listed below. For simplicity, properties are listed as follows:
%
%   f < |x| <- the smooth function f is always less than |x|
%   f >= |x| <- f is always greater than or equal to x
%       ... etc ...
%   f(0) > 0 <- f is greater than 0 at x = 0
%       ... etc ...
%   lim f(x) -> |x| <- f approaches |x| as x -> +/- Inf
%   convex <- f is convex everywhere
%   nonconvex <- f is not convex everywhere
%
% Types include:
%   'sqrt': sqrt(x.^2 + (1/s)^2)
%       Commonly used in optimisation problems (e.g. Kelly 2017)
%       f > |x|     f(0) > 0    lim f(x) -> |x|     convex
%   'tanh': x.*tanh(x*s)
%       Equivalent to Boltzmann operator (max approximation) for two inputs 
%       (-x,x). See also Kelly 2017.
%       f < |x|     f(0) = 0    lim f(x) -> |x|     nonconvex
%   'antilogistic': 1/s*( log(1+exp(-s*x))+log(1+exp(s*x)) )
%       Used by Schmidt2007 (doi:10.1007/978-3-540-74958-5_28) from the 
%       'softplus' function proposed by Chen1996 (doi:10.1007/BF00249052).
%       Based on the antiderivative of the logistic function.
%       f > |x|     f(0) > 0    lim f(x) -> |x|     convex
%   'antilogistic0' 2/s*log(1+exp(s*x)) - x - 2/s*log(2)
%       Similar basis to antilogistic, but offset such that f(0) = 0
%       See https://math.stackexchange.com/q/1115033
%       f < |x|     f(0) = 0    convex

switch lower(type)
    case 'sqrt'
        f_absSmooth = @(x,s) sqrt(x.^2 + (1/s)^2);
    case 'tanh'
        f_absSmooth = @(x,s) x.*tanh(x*s);
    case 'antilogistic'
        f_absSmooth = @(x,s) 1/s*( log(1+exp(-s*x))+log(1+exp(s*x)) );
    case 'antilogistic0'
        f_absSmooth = @(x,s) 2/s*log(1+exp(s*x)) -x - 2/s*log(2);
    otherwise
        error('Type not recognized')
end
        