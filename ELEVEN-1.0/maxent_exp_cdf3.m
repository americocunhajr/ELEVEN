
% -----------------------------------------------------------------
%  maxent_exp_cdf3.m
%
%  This functions computes the CDF function for a MaxEnt
%  truncated exponential distribution with three parameters.
%
%  input:
%  x      - input value
%  lambda - (2 x 1) Lagrange multipliers vector
%  xmin   - support lower bound
%
%  output:
%  cdf_x - CDF function
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Oct 22, 2018
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function cdf_x = maxent_exp_cdf3(x,lambda,xmin)

    % check number of arguments
    if nargin < 3
        error('Too few inputs.')
    elseif nargin > 3
        error('Too many inputs.')
    end
    
    % ensure lambda is a column vector
    lambda = lambda(:);
    
    % check if lambda is a 3-dimensional vector
    if length(lambda) ~= 3
        error('lambda must be a (3 x 1) array')
    end
        
	% Lagrange multipliers
	L0 = lambda(1);
	L1 = lambda(2);
    L2 = lambda(3);
    
    if L2 == 0
        error('lambda 2 must be a nonzero parameter')
    end
    
    % auxiliar variables
    a = exp(-L0 + 0.25*L1^2/L2);
    b = L2;
    c = -0.5*L1/L2;
    
    if b > 0
        
        % auxiliar variables
        A = 0.5*sqrt(pi)*a/sqrt(b);
        xminus = sqrt(b)*(c-xmin);
         xplus = sqrt(b)*(c-x);
        
        % CDF function
        cdf_x = A*(erf(xminus)-erf(xplus));
	
    else
        
        % auxiliar variables
        b = abs(b);
        A = 0.5*sqrt(pi)*a/sqrt(b);
        xminus = sqrt(b)*(c-xmin);
         xplus = sqrt(b)*(c-x);
        
        % CDF function
        cdf_x = A*(erfi(xminus)-erfi(xplus));
    end

return
% -----------------------------------------------------------------
