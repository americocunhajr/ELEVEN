
% -----------------------------------------------------------------
%  elevator_braking_force_robust.m
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
function Fobj = elevator_braking_force_robust(x,param1,param2,MCsmp)

% model parameters
   %Fs = param1(1);
   Fg = param1(2);
   Fb = param1(3);
  mu1 = param1(4);
  mu2 = param1(5);
  mu4 = param1(6);
%alpha = param1(7);
    a = param1(8);
    b = param1(9);
    c = param1(10);
    d = param1(11);
    e = param1(12);
    f = param1(13);
    l = param1(14);
    m = param1(15);
    n = param1(16);
    R = param1(17);

% stochastic solver parameters
Ns = param2(1);
%P_ref = param2(2);
%ylow = param2(3);
%yupp = param2(4);


% design variables
a = x(1);
c = x(2);

% preallocate memory for braking force samples
MC_Fh = zeros(Ns,1);

for imc=1:Ns
    
    % i-th cossine of alpha
    cos_alpha = cos(MCsmp(1,imc));
    
    % i-th cossine of alpha
    sin_alpha = sin(MCsmp(1,imc));
    
    % i-th spring force reaction
    Fs = MCsmp(2,imc);
    
    % i-th normal forces
    N4 = (0.5*(Fg+Fb)*l - Fs*a)/(mu4*(n+l)-m);
    N1 = (N4-(a*mu2*Fs)/(d+e*mu2))/(mu1*sin_alpha+cos_alpha+mu2*(b*mu1-c)/(d+e*mu2));
    N2 = (a*Fs + (b*mu1-c)*N1)/(d+e*mu2);
    N3 = mu2*N2 + (mu1*sin_alpha + cos_alpha)*N1;
    
    % store i-th braking force realization
    MC_Fh(imc,1) = mu1*N1 + mu2*N2 + (f/R)*N3 + mu4*N4;

end

% objective function
beta1 = 0.2;
beta2 = 0.2;
beta3 = 0.2;
beta4 = 0.4;
Fobj = beta1*min(MC_Fh) + beta2*max(MC_Fh) + ...
       beta3*mean(MC_Fh) + beta4*(1/std(MC_Fh));
                          

return
% -----------------------------------------------------------------
