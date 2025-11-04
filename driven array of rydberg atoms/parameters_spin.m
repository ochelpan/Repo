
N=10; % number of spins
ncopies=1;% 
omega0=1.; % level splitting
Omega=1.; % external drive
J0=1; %interaction strength
Alpha=6; %power low decay of interaction with distance 
gamma_d=0.1; %dissipation rate
kappa=.0;% dephasing rate

%% don't touch this parameter!!!
is_noise=1;% 0 when no noise and 1 when noise is added
is_simplified=0; %  we treat spins as fields and total spin is preserved
% is_simplified=1;% we work with simplified eom; eom are linear is spins 

%% initial direction of spin
theta=0;%pi/2;
phi=0;


t_in=0; %initial time
dt=1e-3; %time step
t_fin=20;%finite time
% adjust_mu=50;
adjust_mu=10; % we save data at every 10-th timestep 
% here timestep should be small as the integrating scheme does not preserve
% the total spin and it may increase over time; the smaller is timestep the
% smaller is increase. Shouldn't be a problem if we use spin-preserving
% integration schemes

% adjust_mu=100;

n_iter=20; % iteration within integrator to reach convergence
epsilon_max=1e-12; 
M=round((t_fin-t_in)/dt)+1; %number of time steps
Mev=(t_fin-t_in)/(dt*adjust_mu); %number of time steps for which we save the data

Ntrc=400;%number of trajectories


% some parameters which are likely here for older versions of code
dim=1;
lambda=1;
d=0.2*lambda;
Nx=1;
Ny=1;
