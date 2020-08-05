
% -----------------------------------------------------------------
%  main_elevator_braking_mc.m
%
%  This script is the main file for a program that simulates
%  the stochastic mechanics, using Monte Carlo method,
%  of an elevator braking system.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Jan 29, 2018
% -----------------------------------------------------------------

clc
clear all
close all



% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' ---------------------------------------------------')
disp(' Elevator Breaking Stochastic Mechanics             ')
disp(' (uncertanty propagation calculation)               ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Piotr Wolszczak                                    ')
disp(' p.wolszczak@pollub.pl                              ')
disp('                                                    ')
disp(' Pawel Lonkwic                                      ')
disp(' plonkwic@gmail.com                                 ')
disp('                                                    ')
disp(' Americo Cunha Jr                                   ')
disp(' americo.cunhajr@gmail.com                          ')
disp('                                                    ')
disp(' Grzegorz Litak                                     ')
disp(' g.litak@pollub.pl                                  ')
disp(' ---------------------------------------------------')
disp('                                                    ')
% -----------------------------------------------------------


% simulation information
% -----------------------------------------------------------
case_name = 'elevator_mc';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------



% define deterministic physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining physical parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

Fs    = 50.0e3;  % spring reaction force (N)
alpha = 37.0;    % nominal cam mechanism angle (degree)
% -----------------------------------------------------------


% generate random parameters
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- generating random parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

% RNG seed
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

% number of random samples
Ns = 10e3;

% angle limit values (rad)
alpha_min = 20.0*(pi/180);
alpha_max = 40.0*(pi/180);

% angle mean value (rad)
mu_alpha = alpha*(pi/180);

% angle dispersion parameter (dimensionless)
delta_alpha = 0.1;

% angle std deviation (rad)
sigma_alpha = delta_alpha*mu_alpha;

% number of mesh points for random variable support
Nx = 1000;

% Lagrange multipliers for MaxEnt distribution
[lambda,pdf_alpha,supp_alpha] = ...
    maxent_lagrange_mc(alpha_min,alpha_max,Nx,[1; mu_alpha]);
l0 = lambda(1);
l1 = lambda(2);

% CDF of Alpha
cdf_alpha = @(x) exp(-l0)*(exp(-l1*alpha_min)-exp(-l1*x))/l1;

% quantile function of Alpha
cdf_inv_alpha = @(x) -log(exp(-l1*alpha_min)-l1*x*exp(l0))/l1;

% uniform samples in [0,1]
U = rand(Ns,1);

% MaxEnt random samples for alpha (inverse transform method)
Alpha = cdf_inv_alpha(U);

% Uniform random samples for alpha
%Alpha = alpha_min + (alpha_max-alpha_min)*rand(Ns,1);

toc
% -----------------------------------------------------------



% Monte Carlo simulation
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- Monte Carlo Simulation --- ');
disp(' ');


% preallocate memory for breaking force samples
MC_Fh   = zeros(Ns,1);

% auxiliar constants
a1 = 0.2056*Fs;
a2 = 10142847/392138000;
a5 = 148.0783;

for imc=1:Ns
    
    if mod(imc,round(sqrt(Ns))) == 0
        disp('')
        disp(imc)
    end
    
    % i-th cossine of alpha
    cos_alpha = cos(Alpha(imc));
    
    % i-th cossine of alpha
    sin_alpha = sin(Alpha(imc));
    
    % auxiliar constants
    a3 = (0.3307*Fs-987.1888)/(cos_alpha+0.1*sin_alpha-0.1258);
    a4 = (cos_alpha-0.1*sin_alpha)/2900;
    
    % breaking force
    Fh_imc = a1 - a2*a3 + a4*a3 - a5;
    
    % store i-th breaking force realization
    MC_Fh(imc,1) = Fh_imc;

end

% MC convergence metric (2nd moment of Fh)
MC_conv = cumsum(MC_Fh.^2)./(1:Ns)';

toc
% -----------------------------------------------------------




% compute the statistics
% -----------------------------------------------------------
tic

disp(' ')
disp(' --- computing statistics --- ');
disp(' ');
disp('    ... ');
disp(' ');

% sample average
alpha_smp_avg = mean(Alpha);
   Fh_smp_avg = mean(MC_Fh);
     
% standard deviation
alpha_std = std(Alpha);
   Fh_std = std(MC_Fh);

% skewness
alpha_skew = skewness(Alpha);
   Fh_skew = skewness(MC_Fh);
     
% kurtosis
alpha_kurt = kurtosis(Alpha)-3;
   Fh_kurt = kurtosis(MC_Fh)-3;

% confidence probability (percentual)
Pc = 95;

% upper percentil
r_plus = 0.5*(100 + Pc);

% lower percentil
r_minus = 0.5*(100 - Pc);

% confidence band upper bound
alpha_upp = prctile(Alpha,r_plus);
   Fh_upp = prctile(MC_Fh,r_plus);

% confidence band lower bound
alpha_low = prctile(Alpha,r_minus);
   Fh_low = prctile(MC_Fh,r_minus);
    
% number of bins
Nbins = round(sqrt(Ns));

% number of KSD points
Nksd = 2*Nbins;

%  histogram estimator
[alpha_bins,alpha_freq] = randvar_pdf(Alpha,Nbins);
[   Fh_bins,   Fh_freq] = randvar_pdf(MC_Fh,Nbins);

% kernal density estimator
[alpha_ksd,alpha_supp] = randvar_ksd(Alpha,Nksd);
[   Fh_ksd,   Fh_supp] = randvar_ksd(MC_Fh,Nksd);

% probability of |X| <= x
%X_plus   = ;
%X_minus  = ;

%P_X_plus   = randvar_probval(X_bins,X_freq,X_plus );
%P_X_minus  = randvar_probval(X_bins,X_freq,X_minus);
%P_X        = P_X_plus - P_X_minus;

toc
% -----------------------------------------------------------



% save simulation data into a file
% -----------------------------------------------------------
tic

disp(' ')
disp(' --- Saving Workspace --- ');
disp(' ');
disp('    ... ');
disp(' ');

save([num2str(case_name),'.mat']);

toc
% -----------------------------------------------------------




% post processing
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');

% plot cam mechanism angle statistics
% ...........................................................
gtitle = ' ';
xlab   = ' cam mechanism angle (degree)';
ylab   = ' probability density function';
xmin   = alpha_min*(180/pi);
xmax   = alpha_max*(180/pi);
ymin   = 0.0;
ymax   = max(1.2*alpha_ksd);
leg1   = 'PDF';
leg2   = 'mean';
leg3   = 'mean \pm std';
leg4   = [' ',num2str(Pc),'% prob.'];
gname  = 'alpha_pdf_fig';
flag   = 'eps';
fig1 = graph_randvar_statistics(supp_alpha*(180/pi),...
                                pdf_alpha,...
                                alpha_smp_avg*(180/pi),...
                                alpha_std*(180/pi),...
                                alpha_low*(180/pi),...
                                alpha_upp*(180/pi),...
                                gtitle,xlab,ylab,...
                                leg1,leg2,leg3,leg4,...
                                xmin,xmax,ymin,ymax,gname,flag);
%close(fig1);
% ...........................................................


% plot breaking force statistics
% ...........................................................
gtitle = ' ';
xlab   = ' braking force (N)';
ylab   = ' probability density function';
xmin   = min(MC_Fh);
xmax   = max(MC_Fh);
ymin   = 0.0;
ymax   = max(1.2*Fh_ksd);
leg1   = 'PDF';
leg2   = 'mean';
leg3   = 'mean \pm std';
leg4   = [' ',num2str(Pc),'% prob.'];
gname  = 'Fh_pdf_fig';
flag   = 'eps';
fig2 = graph_randvar_statistics(Fh_supp,Fh_ksd,Fh_smp_avg,...
                                Fh_std,Fh_low,Fh_upp,...
                                gtitle,xlab,ylab,...
                                leg1,leg2,leg3,leg4,...
                                xmin,xmax,ymin,ymax,gname,flag);
%close(fig2);
% ...........................................................


% plot cam mechanism angle samples
% ...........................................................
gtitle = ' ';
xlab   = ' sample index';
ylab   = ' cam mechanism angle (degree)';
xmin   = 0;
xmax   = Ns;
ymin   = 0.95*alpha_min*(180/pi);
ymax   = 1.05*alpha_max*(180/pi);
leg1   = 'sample';
leg2   = 'mean';
leg3   = 'mean \pm std';
leg4   = [' ',num2str(Pc),'% prob.'];
gname  = 'alpha_samples_fig';
flag   = 'eps';
fig3 = graph_randvar_samples([1:Ns],Alpha*(180/pi),...
                             alpha_smp_avg*(180/pi),...
                             alpha_std*(180/pi),...
                             alpha_low*(180/pi),...
                             alpha_upp*(180/pi),...
                             gtitle,xlab,ylab,...
                             leg1,leg2,leg3,leg4,...
                             xmin,xmax,ymin,ymax,gname,flag);
%close(fig3);
% ...........................................................


% plot breaking force samples
% ...........................................................
gtitle = ' ';
xlab   = ' sample index';
ylab   = ' braking force (N)';
xmin   = 0;
xmax   = Ns;
ymin   = 0.999*min(MC_Fh);
ymax   = 1.001*max(MC_Fh);
leg1   = 'sample';
leg2   = 'mean';
leg3   = 'mean \pm std';
leg4   = [' ',num2str(Pc),'% prob.'];
gname  = 'Fh_samples_fig';
flag   = 'eps';
fig4 = graph_randvar_samples([1:Ns],MC_Fh,Fh_smp_avg,...
                             Fh_std,Fh_low,Fh_upp,...
                             gtitle,xlab,ylab,...
                             leg1,leg2,leg3,leg4,...
                             xmin,xmax,ymin,ymax,gname,flag);
%close(fig4);
% ...........................................................


% plot MC convergence metric
% ...........................................................
gtitle = ' ';
xlab   = ' number of MC realizations';
ylab   = ' convergence metric';
xmin   = 0;
xmax   = Ns;
ymin   = 'auto';
ymax   = 'auto';
gname  = 'MC_conv';
flag   = 'eps';

fig5   = graph_type1((1:Ns),MC_conv,gtitle,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig5);
% ...........................................................
                                    
toc
% -----------------------------------------------------------


