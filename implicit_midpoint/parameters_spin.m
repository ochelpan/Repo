% folder to save data for each trajectory
foldername_date=strcat('/Users/oksanachelpanova/Matlab_code/parallel_code/');
[status, msg, msgID] = mkdir(foldername_date); % create folder to save data if it does not exist

DO_PARALLEL=0; % set to 1 if tun on GPU
ncopies=10; %doing calculations in bunches 
Nrep=500;
Ntrc=ncopies*Nrep;%number of realizations
Knoise=5000;% sampling noise for the next Knoise steps at once


N=1; % number of spins
omega0=0.; % level splitting
Omega=0.0;% coherent drive

Eps=(1e-3); % perturbation of the steady state to calculate two-time correlation function

%% don't touch this parameter!!!
is_noise=1;% 0 when no noise and 1 when noise is added
is_simplified=0; %  we treat spins as fields and total spin is preserved
% is_simplified=1;% we work with simplified eom; eom are linear is spins 

%% initial direction of spin
% theta=pi/2;
% phi=0;
theta=0.;%1e-2;%pi/2;
phi=0.;%1e-1;

tau0=1e-1;% duration of the steady state preparation 
t_in=0.; %initial time

 t_fin=50.;%finite time 

dt=1e-1;%timestep 
adjust_mu=1;% i save every 10-th value of spins (every dt_eff=0.01)

n_iter=8;
epsilon_max=1e-10;
M=round((t_fin-t_in)/dt)+1; %number of time steps
Mev=(t_fin-t_in)/(dt*adjust_mu); %number of time steps that we save
Mev0=tau0/(dt*adjust_mu);
M0=round((tau0)/dt)+1;%number of timesteps at the preparation stage


% Ntrc=floor(50e3/N);%number of realizations


ind_check=10;%floor(0.1*Ntrc/N);




dim=1;
lambda=1;

d=0.2*lambda;
Nx=1;
Ny=1;
J0=0;
% Alpha=6;

gamma_d=1.; %dissipation rate
gamma_u=1.;
kappa=.0;% dephasing rate