function [t,S_m,S_z] = solver_spins_Rydbergs_real_loop(sigma_m_0,sigma_z_0,J,h)

% n_iter=20;
% epsilon_max=1e-6;
parameters_spin;
n=1;
mu=adjust_mu-1; % here, spin value are saved for every adjust_mu time step 
L=length(sigma_m_0);

%initializing output array
S_m(1:L,1:Mev)=0;
S_z(1:L,1:Mev)=0;
t(1:Mev)=0;

%% preparation
rng('shuffle','v5normal');

% sampling noise; In fact can also be done directly
% within the loop, which saves a lot of memory

noise_gamma_diag1=is_noise*(randn(L,M))*sqrt(gamma_d*dt); % noise in diagonal basis
noise_gamma_diag2=is_noise*(randn(L,M))*sqrt(gamma_d*dt); % noise in diagonal basis
noise_kappa_diag=is_noise*(randn(L,M))*sqrt(4*kappa*dt); % noise in diagonal basis

noise_gamma_real=noise_gamma_diag1;
noise_kappa_real=noise_kappa_diag;
noise_gamma_imag=noise_gamma_diag2;

%initial conditions

sx(1:L,1)=2*real(sigma_m_0);
sy(1:L,1)=2*imag(sigma_m_0);
sz(1:L,1)=sigma_z_0;

% evolution
for m=2:M
    mu=mu+1;
    % noise at each time step
    noiseGammaPrev_real = noise_gamma_real(:, m - 1); 
    noiseKappaPrev_real = noise_kappa_real(:, m - 1);    
    noiseGammaPrev_imag = noise_gamma_imag(:, m - 1);
 %% previous step

    A11=-2*h.*sy-0.5*sy.*(J*sz)+gamma_d/2*sx.*sz;
    A12=-2*Omega*sz+2*h.*sx+0.5*sx.*(J*sz)+gamma_d/2*sy.*sz;
    A13=2*Omega*sy-gamma_d/2*(sx.*sx+sy.*sy);
     %% new step predictor
    sx_tilde=sx+A11*dt+noiseGammaPrev_imag.*sz+noiseKappaPrev_real.*sy;
    sy_tilde=sy+A12*dt+noiseGammaPrev_real.*sz-noiseKappaPrev_real.*sx;
    sz_tilde=sz+A13*dt-(noiseGammaPrev_real.*sy+noiseGammaPrev_imag.*sx); 

    for i=1:n_iter
    
    A11_tilde=-2*h.*sy_tilde-0.5*sy_tilde.*(J*sz_tilde)+gamma_d/2*sx_tilde.*sz_tilde;
    A12_tilde=-2*Omega*sz_tilde+2*h.*sx_tilde+0.5*sx_tilde.*(J*sz_tilde)+gamma_d/2*sy_tilde.*sz_tilde;
    A13_tilde=2*Omega*sy_tilde-gamma_d/2*(sx_tilde.*sx_tilde+sy_tilde.*sy_tilde);


    sx_new=sx+0.5*(A11+A11_tilde)*dt+0.5*noiseGammaPrev_imag.*(sz+sz_tilde)+0.5*noiseKappaPrev_real.*(sy+sy_tilde);
    sy_new=sy+0.5*(A12+A12_tilde)*dt+0.5*noiseGammaPrev_real.*(sz+sz_tilde)-0.5*noiseKappaPrev_real.*(sx+sx_tilde);
    sz_new=sz+0.5*(A13+A13_tilde)*dt-0.5*(noiseGammaPrev_real.*(sy+sy_tilde)+noiseGammaPrev_imag.*(sx+sx_tilde)); 


    % breaking loop if precision is high
    if max(abs(sz_new-sz_tilde))<epsilon_max
        if max(abs(sy_new-sy_tilde))<epsilon_max
            if max(abs(sx_new-sx_tilde))<epsilon_max
%                 i
                break;
            end
        end
    end

    sx_tilde=sx_new;
    sy_tilde=sy_new;
    sz_tilde=sz_new;

    end

    % saving data 
    if mu==adjust_mu
       S_m(1:L,n)= 0.5*(sx_new-1i*sy_new);
       S_z(1:L,n)= sz_new;
       t(n)=t_in+(m-1)*dt;

        mu=0;
        n=n+1;
    end
    
    sx=sx_new;
    sy=sy_new;
    sz=sz_new;

end

end