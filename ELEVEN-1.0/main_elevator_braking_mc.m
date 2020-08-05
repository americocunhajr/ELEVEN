
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
%  last update: Nov 4, 2018
% -----------------------------------------------------------------

clc
clear
close all



% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' ---------------------------------------------------')
disp(' Elevator Braking Stochastic Mechanics              ')
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



% physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining physical parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

Fs    = 42.0e3;  % spring reaction force (N)
Fg    = 6.0e3;  % lift cabin weight (N)
Fb    = 1.2e3;  % lift cabin inertial force (N)

alpha = 6.0;    % nominal cam mechanism angle (degree)

a = 55.00e-3;  % geometric dimension (m)
b = 16.59e-3;  % geometric dimension (m)
c = 52.70e-3;  % geometric dimension (m)
d = 34.50e-3;  % geometric dimension (m)
e = 60.66e-3;  % geometric dimension (m)
f = 0.005e-3;  % geometric dimension (m)
l = 49.00e-3;  % geometric dimension (m)
m = 40.00e-3;  % geometric dimension (m)
n = 17.50e-3;  % geometric dimension (m)
R = 29.00e-3;  % geometric dimension (m)

mu1 = 0.10;   % friction coeficient
mu2 = 0.10;   % friction coeficient
mu4 = 0.15;   % friction coeficient
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
Ns = 4096;

% angle limit values (rad)
alpha_min =  0.0*(pi/180);
alpha_max = 18.0*(pi/180);

% angle mean value (rad)
mu_alpha = 6*(pi/180);


% Fs limit values (N)
Fs_min = 0.0;
Fs_max = 56.0e3;

% Fs mean value (N)
mu_Fs = 42.0e3;

% number of mesh points for random variable support
Nx = 1000;

% Lagrange multipliers for MaxEnt distribution
[lambda_alpha,pdf_alpha,supp_alpha] = ...
    maxent_lagrange_mc(alpha_min,alpha_max,Nx,[1; mu_alpha]);
[lambda_Fs,pdf_Fs,supp_Fs] = ...
    maxent_lagrange_mc(Fs_min,Fs_max,Nx,[1; mu_Fs]);

% CDF for MaxEnt distribution
cdf_alpha = maxent_exp_cdf2(supp_alpha,lambda_alpha,alpha_min);
cdf_Fs    = maxent_exp_cdf2(supp_Fs   ,lambda_Fs   ,Fs_min);

% uniform samples in [0,1]
U = rand(Ns,1);

% random samples for alpha (inverse transform method)
MC_alpha = maxent_exp_quant2(U(:,1),lambda_alpha,alpha_min);
MC_Fs    = maxent_exp_quant2(U(:,1),lambda_Fs   ,Fs_min);

toc
% -----------------------------------------------------------



% Monte Carlo simulation
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- Monte Carlo Simulation --- ');
disp(' ');


% preallocate memory for braking force samples
MC_Fh = zeros(Ns,1);

for imc=1:Ns
    
    if mod(imc,round(sqrt(Ns))) == 0
        disp('')
        disp(imc)
    end
    
    % i-th cossine of alpha
    cos_alpha = cos(MC_alpha(imc));
    
    % i-th cossine of alpha
    sin_alpha = sin(MC_alpha(imc));
    
    % i-th spring force reaction
    Fs = MC_Fs(imc);
    
    % i-th normal forces
    N4 = (0.5*(Fg+Fb)*l - Fs*a)/(mu4*(n+l)-m);
    N1 = (N4-(a*mu2*Fs)/(d+e*mu2))/(mu1*sin_alpha+cos_alpha+mu2*(b*mu1-c)/(d+e*mu2));
    N2 = (a*Fs + (b*mu1-c)*N1)/(d+e*mu2);
    N3 = mu2*N2 + (mu1*sin_alpha + cos_alpha)*N1;
    
    % store i-th braking force realization
    MC_Fh(imc,1) = mu1*N1 + mu2*N2 + (f/R)*N3 + mu4*N4;

end

% MC convergence metrics (low order statistics)
MC_conv_Fh_mean = cumsum(MC_Fh)./(1:Ns)';
MC_conv_Fh_2nd  = cumsum(MC_Fh.^2)./(1:Ns)';
MC_conv_Fh_std  = sqrt(MC_conv_Fh_2nd - MC_conv_Fh_mean.^2);


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
alpha_mean = mean(MC_alpha);
   Fs_mean = mean(MC_Fs);
   Fh_mean = mean(MC_Fh);
     
% standard deviation
alpha_std = std(MC_alpha);
   Fs_std = std(MC_Fs);
   Fh_std = std(MC_Fh);

% skewness
alpha_skew = skewness(MC_alpha);
   Fs_skew = skewness(MC_Fs);
   Fh_skew = skewness(MC_Fh);
     
% kurtosis
alpha_kurt = kurtosis(MC_alpha)-3;
   Fs_kurt = kurtosis(MC_Fs)-3;
   Fh_kurt = kurtosis(MC_Fh)-3;

% confidence probability (percentual)
Pc = 95;

% upper percentil
r_plus = 0.5*(100 + Pc);

% lower percentil
r_minus = 0.5*(100 - Pc);

% confidence band upper bound
alpha_upp = prctile(MC_alpha,r_plus);
   Fs_upp = prctile(MC_Fs,r_plus);
   Fh_upp = prctile(MC_Fh,r_plus);

% confidence band lower bound
alpha_low = prctile(MC_alpha,r_minus);
   Fs_low = prctile(MC_Fs,r_minus);
   Fh_low = prctile(MC_Fh,r_minus);
    
% number of bins
Nbins = round(sqrt(Ns));

% number of KSD points
Nksd = 2*Nbins;

%  histogram estimator
[bins_alpha,freq_alpha] = randvar_pdf(MC_alpha,Nbins);
[bins_Fs   ,freq_Fs   ] = randvar_pdf(MC_Fs   ,Nbins);
[bins_Fh   ,freq_Fh   ] = randvar_pdf(MC_Fh   ,Nbins);

% kernal density estimator
min_Fh = floor(min(MC_Fh));
max_Fh = ceil(max(MC_Fh));
%[ksd_Fh,supp_Fh] = randvar_ksd(MC_Fh,Nksd);
[ksd_Fh,supp_Fh] = ksdensity(MC_Fh,bins_Fh,'Support',[min_Fh max_Fh]);
 ksd_Fh = ksd_Fh(2:end-1);
supp_Fh = supp_Fh(2:end-1);

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
xlab   = ' cam angle (degree)';
ylab   = ' probability density function';
xmin   = alpha_min*(180/pi);
xmax   = alpha_max*(180/pi);
ymin   = 0.0;
ymax   = max(1.2*pdf_alpha);
leg1   = 'PDF';
leg2   = 'mean';
leg3   = 'mean \pm std';
leg4   = [' ',num2str(Pc),'% prob.'];
gname  = 'alpha_pdf_fig';
flag   = 'eps';
fig1 = graph_randvar_statistics(supp_alpha*(180/pi),...
                                pdf_alpha,...
                                alpha_mean*(180/pi),...
                                alpha_std*(180/pi),...
                                alpha_low*(180/pi),...
                                alpha_upp*(180/pi),...
                                bins_alpha*(180/pi),freq_alpha,...
                                gtitle,xlab,ylab,...
                                leg1,leg2,leg3,leg4,...
                                xmin,xmax,ymin,ymax,gname,flag);
%close(fig1);
% ...........................................................



% plot spring force statistics
% ...........................................................
gtitle = ' ';
xlab   = ' spring force (N)';
ylab   = ' probability density function';
xmin   = min(MC_Fs);
xmax   = max(MC_Fs);
ymin   = 0.0;
ymax   = max(1.2*pdf_Fs);
leg1   = 'PDF';
leg2   = 'mean';
leg3   = 'mean \pm std';
leg4   = [' ',num2str(Pc),'% prob.'];
gname  = 'Fs_pdf_fig';
flag   = 'eps';
fig2 = graph_randvar_statistics(supp_Fs,pdf_Fs,Fs_mean,...
                                Fs_std,Fs_low,Fs_upp,...
                                bins_Fs,freq_Fs,...
                                gtitle,xlab,ylab,...
                                leg1,leg2,leg3,leg4,...
                                xmin,xmax,ymin,ymax,gname,flag);
%close(fig2);
% ...........................................................


% plot braking force statistics
% ...........................................................
gtitle = ' ';
xlab   = ' braking force (N)';
ylab   = ' probability density function';
%xmin   = Fh_low;
%xmax   = Fh_upp;
xmin   = min(MC_Fh);
xmax   = max(MC_Fh);
ymin   = 0.0;
ymax   = max(1.2*ksd_Fh);
leg1   = 'PDF';
leg2   = 'mean';
leg3   = 'mean \pm std';
leg4   = [' ',num2str(Pc),'% prob.'];
gname  = 'Fh_pdf_fig';
flag   = 'eps';
fig3 = graph_randvar_statistics(supp_Fh,ksd_Fh,Fh_mean,...
                                Fh_std,Fh_low,Fh_upp,...
                                bins_Fh,freq_Fh,...
                                gtitle,xlab,ylab,...
                                leg1,leg2,leg3,leg4,...
                                xmin,xmax,ymin,ymax,gname,flag);
%close(fig3);
% ...........................................................


% plot cam mechanism angle samples
% ...........................................................
gtitle = ' ';
xlab   = ' sample index';
ylab   = ' cam angle (degree)';
xmin   = 0;
xmax   = Ns;
ymin   = min(MC_alpha)*(180/pi);
ymax   = max(MC_alpha)*(180/pi);
leg1   = 'sample';
leg2   = 'mean';
leg3   = 'mean \pm std';
leg4   = [' ',num2str(Pc),'% prob.'];
gname  = 'alpha_samples_fig';
flag   = 'eps';
fig4 = graph_randvar_samples([1:Ns],MC_alpha*(180/pi),...
                             alpha_mean*(180/pi),...
                             alpha_std*(180/pi),...
                             alpha_low*(180/pi),...
                             alpha_upp*(180/pi),...
                             gtitle,xlab,ylab,...
                             leg1,leg2,leg3,leg4,...
                             xmin,xmax,ymin,ymax,gname,flag);
%close(fig4);
% ...........................................................


% plot spring force samples
% ...........................................................
gtitle = ' ';
xlab   = ' sample index';
ylab   = ' spring force (N)';
xmin   = 0;
xmax   = Ns;
ymin   = min(MC_Fs);
ymax   = max(MC_Fs);
leg1   = 'sample';
leg2   = 'mean';
leg3   = 'mean \pm std';
leg4   = [' ',num2str(Pc),'% prob.'];
gname  = 'Fs_samples_fig';
flag   = 'eps';
fig5 = graph_randvar_samples([1:Ns],MC_Fs,Fs_mean,...
                             Fs_std,Fs_low,Fs_upp,...
                             gtitle,xlab,ylab,...
                             leg1,leg2,leg3,leg4,...
                             xmin,xmax,ymin,ymax,gname,flag);
%close(fig5);
% ...........................................................



% plot breaking force samples
% ...........................................................
gtitle = ' ';
xlab   = ' sample index';
ylab   = ' braking force (N)';
xmin   = 0;
xmax   = Ns;
ymin   = min(MC_Fh);
ymax   = max(MC_Fh);
leg1   = 'sample';
leg2   = 'mean';
leg3   = 'mean \pm std';
leg4   = [' ',num2str(Pc),'% prob.'];
gname  = 'Fh_samples_fig';
flag   = 'eps';
fig6 = graph_randvar_samples([1:Ns],MC_Fh,Fh_mean,...
                             Fh_std,Fh_low,Fh_upp,...
                             gtitle,xlab,ylab,...
                             leg1,leg2,leg3,leg4,...
                             xmin,xmax,ymin,ymax,gname,flag);
%close(fig6);
% ...........................................................


% plot MC convergence Fh std
% ...........................................................
gtitle = ' ';
xlab   = ' number of samples';
ylab   = ' braking force std deviation';
xmin   = 0;
xmax   = Ns;
ymin   = 'auto';
ymax   = 'auto';
gname  = 'MC_conv_Fh_std';
flag   = 'eps';

fig7   = graph_type1((1:Ns),MC_conv_Fh_std,gtitle,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig7);
% ...........................................................


% plot MC convergence Fh std
% ...........................................................
gtitle = ' ';
xlab   = ' number of samples';
ylab   = ' braking force mean';
xmin   = 0;
xmax   = Ns;
ymin   = 'auto';
ymax   = 'auto';
gname  = 'MC_conv_Fh_mean';
flag   = 'eps';

fig8   = graph_type1((1:Ns),MC_conv_Fh_mean,gtitle,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig8);
% ...........................................................
                                    
toc
% -----------------------------------------------------------


