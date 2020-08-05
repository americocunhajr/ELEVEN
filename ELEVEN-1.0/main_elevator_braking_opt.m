
% -----------------------------------------------------------------
%  main_elevator_braking_opt.m
%
%  This script is the main file for a program that optimizes
%  the mechanics of an elevator braking system.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 2, 2018
% -----------------------------------------------------------------

clc
clear
close all



% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' ---------------------------------------------------')
disp(' Elevator Braking Optimization                      ')
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
case_name = 'elevator_opt';

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

Fs    = 43.0e3;  % spring reaction force (N)
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
U = rand(Ns,2);

% samples for Monte Carlo simulation
MCsmp = zeros(2,Ns);

% random samples for alpha (inverse transform method)
MCsmp(1,:) = maxent_exp_quant2(U(:,1),lambda_alpha,alpha_min);
%MCsmp(1,:) = alpha_min + (alpha_max-alpha_min)*rand(Ns,1);
MCsmp(2,:) = maxent_exp_quant2(U(:,2),lambda_Fs   ,Fs_min);

toc
% -----------------------------------------------------------


% solve optimization problem
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- Optimization --- ');
disp(' ');

% solver optional parameters
opt = optimoptions('fmincon','Display','iter',...
                             'Algorithm','sqp',... 
                             'SpecifyObjectiveGradient' ,false,...
                             'SpecifyConstraintGradient',false);

% reference probability
P_ref = 0.05;

% bounds for Fh
ylow = 0.5e3;
yupp = 5.0e3;
                         
% model parameters vector
param1 = [Fs Fg Fb mu1 mu2 mu4 alpha a b c d e f l m n R];

% stochastic solver parameters vector
param2 = [Ns P_ref ylow yupp];

% design variables limits
min_x1 = 50.0e-3;
max_x1 = 60.0e-3;

min_x2 = 50.0e-3;
max_x2 = 55.0e-3;
                         
% lower and upper bounds
lb = [min_x1 min_x2];
ub = [max_x1 max_x2];

% initial guess
x0 = 0.5*(lb + ub);

% objetive and constraint functions
  ObjF1 = @(x) (-1)*elevator_braking_force(x,param1);
  ObjF2 = @(x) (-1)*elevator_braking_force_robust(x,param1,param2,MCsmp);
ConstrF = @(x) elevator_braking_force_prob(x,param1,param2,MCsmp);

% SQP solver with nonlinear constraint
[x_opt1,Fh_opt1,flag1] = fmincon(ObjF1,x0,[],[],[],[],lb,ub,[],opt);
[x_opt2,Fh_opt2,flag2] = fmincon(ObjF2,x0,[],[],[],[],lb,ub,ConstrF,opt);

% maximum value
Fh_opt1 = (-1)*Fh_opt1;
Fh_opt2 = (-1)*Fh_opt2;

% discretize objective functions (for vizualization)
Nx1 = 30;
Nx2 = 30;

grid_x1 = linspace(min_x1,max_x1,Nx1);
grid_x2 = linspace(min_x2,max_x2,Nx2);
grid_F1 = zeros(Nx1,Nx2);
grid_F2 = zeros(Nx1,Nx2);
grid_G1 = zeros(Nx1,Nx2);
grid_G2 = zeros(Nx1,Nx2);

for nx1=1:Nx1
    for nx2 = 1:Nx2
        
        x = [grid_x1(nx1) grid_x2(nx2)];
        
        [c_ineq,c_eq] = elevator_braking_force_prob(x,param1,param2,MCsmp);

        grid_F1(nx1,nx2) = elevator_braking_force(x,param1);
        grid_F2(nx1,nx2) = elevator_braking_force_robust(x,param1,param2,MCsmp);
        grid_G1(nx1,nx2) = c_ineq(1);
        %grid_G2(nx1,nx2) = c_ineq(2);
    end
end

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



% plot deterministic countour map
% -----------------------------------------------------------
gtitle = ' classical optimization';
xlab   = ' geometric dimension a (mm)';
ylab   = ' geometric dimension c (mm)';
xmin   = min_x1*1.0e3;
xmax   = max_x1*1.0e3;
ymin   = min_x2*1.0e3;
ymax   = max_x2*1.0e3;
gname  = 'elevator_opt_Fh_map';
flag   = 'eps';
fig1   = graph_contourf_pnt(grid_x1*1.0e3,grid_x2*1.0e3,grid_F1',...
                            x_opt1(1)*1.0e3,x_opt1(2)*1.0e3,gtitle,...
                            xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1);
% -----------------------------------------------------------


% plot robust constraint
% -----------------------------------------------------------
gtitle = ' probabilistic constraint';
xlab   = ' geometric dimension a (mm)';
ylab   = ' geometric dimension c (mm)';
xmin   = min_x1*1.0e3;
xmax   = max_x1*1.0e3;
ymin   = min_x2*1.0e3;
ymax   = max_x2*1.0e3;
gname1  = 'elevator_robustopt_Fh_constr1';
gname2  = 'elevator_robustopt_Fh_constr2';
flag   = 'eps';
fig21  = graph_contourf_pnt(grid_x1*1.0e3,grid_x2*1.0e3,grid_G1',...
                            x_opt2(1)*1.0e3,x_opt2(2)*1.0e3,gtitle,...
                            xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
                        
% fig22  = graph_contourf_pnt(grid_x1*1.0e3,grid_x2*1.0e3,grid_G2',...
%                             x_opt2(1)*1.0e3,x_opt2(2)*1.0e3,gtitle,...
%                             xlab,ylab,xmin,xmax,ymin,ymax,gname2,flag);

%close(fig2);
% -----------------------------------------------------------


% plot robust countour map
% -----------------------------------------------------------
gtitle = ' robust optimization';
xlab   = ' geometric dimension a (mm)';
ylab   = ' geometric dimension c (mm)';
xmin   = min_x1*1.0e3;
xmax   = max_x1*1.0e3;
ymin   = min_x2*1.0e3;
ymax   = max_x2*1.0e3;
gname  = 'elevator_robustopt_Fh_map';
flag   = 'eps';
fig3   = graph_contourf_pnt(grid_x1*1.0e3,grid_x2*1.0e3,grid_F2',...
                            x_opt2(1)*1.0e3,x_opt2(2)*1.0e3,gtitle,...
                            xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig3);
% -----------------------------------------------------------

                                 
toc
% -----------------------------------------------------------


