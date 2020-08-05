
% -----------------------------------------------------------------
%  maxent_exp_quant2.m
%
%  This functions computes the quantile function for a MaxEnt
%  truncated exponential distribution with two parameters.
%
%  input:
%  x      - input value
%  lambda - (2 x 1) Lagrange multipliers vector
%  xmin   - support lower bound
%
%  output:
%  Q  - quantile function
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Oct 18, 2018
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function Q = maxent_exp_quant2(x,lambda,xmin)

    % check number of arguments
    if nargin < 3
        error('Too few inputs.')
    elseif nargin > 3
        error('Too many inputs.')
    end
    
    % ensure lambda is a column vector
    lambda = lambda(:);
    
    % check if lambda is a 2-dimensional vector
    if length(lambda) ~= 2
        error('lambda must be a (2 x 1) array')
    end
    
	% Lagrange multipliers
	L0 = lambda(1);
	L1 = lambda(2);
    
	% quantile function
    if L1 ~= 0
        Q = -log(abs(exp(- L1*xmin) - L1*exp(L0)*x))/L1;
    else
        Q = x*exp(L0) + xmin;
    end

return
% -----------------------------------------------------------------
