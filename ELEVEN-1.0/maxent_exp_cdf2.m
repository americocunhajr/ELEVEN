
% -----------------------------------------------------------------
%  maxent_exp_cdf2.m
%
%  This functions computes the CDF function for a MaxEnt
%  truncated exponential distribution with two parameters.
%
%  input:
%  supp_x - support vector or value
%  lambda - (2 x 1) Lagrange multipliers vector
%  xmin   - support lower bound
%
%  output:
%  cdf_x - CDF function
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Oct 18, 2018
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function cdf_x = maxent_exp_cdf2(x,lambda,xmin)

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
    
	% CDF function
    if L1 ~= 0
        cdf_x = exp(-L0)*(exp(-L1*xmin)-exp(-L1*x))/L1;
    else
        cdf_x = exp(-L0)*(x-xmin);
    end

return
% -----------------------------------------------------------------
