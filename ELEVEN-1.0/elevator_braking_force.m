
% -----------------------------------------------------------------
%  elevator_braking_force.m
%
%  This function defines the objective function for
%  fmincon routine.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 2, 2018
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function Fh = elevator_braking_force(x,param)

% model parameters
   Fs = param(1);
   Fg = param(2);
   Fb = param(3);
  mu1 = param(4);
  mu2 = param(5);
  mu4 = param(6);
alpha = param(7);
    a = param(8);
    b = param(9);
    c = param(10);
    d = param(11);
    e = param(12);
    f = param(13);
    l = param(14);
    m = param(15);
    n = param(16);
    R = param(17);

% design variables
a = x(1);
c = x(2);

% auxiliar variables
sin_alpha = sin(alpha);
cos_alpha = cos(alpha);

% normal forces
N4 = (0.5*(Fg+Fb)*l - Fs*a)/(mu4*(n+l)-m);
N1 = (N4-(a*mu2*Fs)/(d+e*mu2))/(mu1*sin_alpha+cos_alpha+mu2*(b*mu1-c)/(d+e*mu2));
N2 = (a*Fs + (b*mu1-c)*N1)/(d+e*mu2);
N3 = mu2*N2 + (mu1*sin_alpha + cos_alpha)*N1;
    
% braking force
Fh = mu1*N1 + mu2*N2 + (f/R)*N3 + mu4*N4;

return
% -----------------------------------------------------------------
