
% -----------------------------------------------------------------
%  maxent_lagrange_mc.m
%
%  This function computes the Lagrange multipliers for a
%  maximum entropy (MaxEnt) distribution given the first
%  N statistical moments constraints.
%
%  Reference:
%  A Matlab Program to Calculate the Maximum Entropy Distributions
%  Author: Ali Mohammad-Djafari
%  In: Smith C.R., Erickson G.J., Neudorfer P.O. (eds)
%  Maximum Entropy and Bayesian Methods, pp 221-233
%  Springer, 1992
%
%  input:
%  xmin    - random variable supp_xort left extreme
%  xmin    - random variable supp_xort right extreme
%  Nx      - number of points for supp_xort discretization
%  mu      - (N x 1) statistical moments
%  lambda0 - (N x 1) Lagrange multipliers initial guess (optional)
%  tol     - numerical tolerance (optional)
%  maxiter - maximum number of iterations (optional)
%
%  output:
%  lambda_x  - (N  x 1) Lagrange multipliers vector
%  pdf_x     - (Nx x 1) random variable PDF
%  entropy_x - (N  x 1) entropy value of the PDF
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Oct 22, 2018
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [lambda_x,pdf_x,supp_x,entropy_x] = ...
             maxent_lagrange_mc(xmin,xmax,Nx,mu,lambda0,tol,maxiter)

    % check number of arguments
    if nargin < 4
        error('Too few inputs.')
    elseif nargin > 7
        error('Too many inputs.')
    end
    
    % check arguments
    if xmin >= xmax
        error('xmin must be less than xmax.')
    end
    
    if Nx < 2
        error('Nx must be an integer greater than one.')
    end
    
    % number of constraints
    N = length(mu);
    
    if nargin == 4
        % prealocate memory for lambda0
        lambda0 = zeros(N,1);
        % initial guess for lambda
        lambda0(1) = log(xmax-xmin);
        % tolerance
        tol = 1.0e-6;
        % maximum of iteration
        maxiter = 20;
    elseif nargin == 5
        % tolerance
        tol = 1.0e-6;
        % maximum of iteration
        maxiter = 20;
    elseif nargin == 6
        % maximum of iteration
        maxiter = 20;
    end
    
    % check arguments
    if length(mu) ~= length(lambda0)
        error('mu and lambda0 vectors must be same length')
    end
    
    % prealocate memory for entropy
    entropy_x = zeros(maxiter,1);
    
    % prealocate memory for phi
	phi = zeros(Nx,2*N-1);
    
    % prealocate memory for G
    G = zeros(2*N-1,1);
    
    % prealocate memory for gnk
    gnk = zeros(N,N);
    
    % discretization of random variable support
    supp_x = linspace(xmin,xmax,Nx)';
    dx   = supp_x(2) - supp_x(1);
    
    % define constraint functions
    % (phi_0(x) = 1 and phi_n(x) = x.^n, n = 1...N)
    phi(:,1) = ones(Nx,1);
    for n = 2:2*N-1
        phi(:,n) = supp_x.*phi(:,n-1);
    end

    % initialize iteration counter
    iter = 0;
    
    % initial guess for lambda
    lambda_x = lambda0;

    % Newton method iteration
    while iter < maxiter
        
        % update iteration counter
        iter = iter + 1;
        
        % compute MaxEnt PDF
        pdf_x = exp(-(phi(:,1:N)*lambda_x));
        
        % compute nonlinear equations G_n(lambda) = mu_n
        for n = 1:2*N-1
            G(n) = dx*sum(phi(:,n).*pdf_x);
        end
        
        % compute entropy value
        entropy_x(iter) = lambda_x'*G(1:N);
        
        % compute Hankel matrix
        for i = 1:N
            gnk(:,i) = -G(i:N+i-1);
        end

        % compute v
        v = mu - G(1:N);
        
        % compute delta
        %delta = gnk\v;
        
        [Q_gnk,R_gnk] = qr(gnk); % QR decomp due to ill-cond of gnk
        delta = R_gnk\(Q_gnk')*v;
        
        %[U,S,V] = svd(gnk);
        %delta = V*((U'*v)./diag(S));
        
        % compute lambda
        lambda_x = lambda_x + delta;
        
        % check convergence
        if abs(delta./lambda_x) < tol
            break
        end
        
        % check convergence
        if iter > 2
            if abs((entropy_x(iter)-entropy_x(iter-1))/entropy_x(iter)) < tol
                break
            end
        end
        
    end
    
    % eliminate numerical resiuduals
    lambda_x = lambda_x.*(abs(lambda_x) > 100*eps);

    % compute PDF
    pdf_x = exp(-(phi(:,1:N)*lambda_x));
    
    % compute the area
    A = trapz(supp_x,pdf_x);
    
    % normalize the PDF
    pdf_x = pdf_x/A;
    
    % correction for lambda_0
    lambda_x(1) = lambda_x(1) + log(1/A);
    
end
