function [S_m,S_z]= solver_spins_corr_emission_real_loop_prep(sigma_m_0,sigma_z_0,J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p)
% this code does stochastic evolution and returns finale state at time
% tau_0 . In principle, one can also save intermediate values of S_m, S_z,
% also t as it is done in solver_spins_corr_emission_real_loop_par(). 

parameters_spin;
L=length(sigma_m_0);


%% preparation
rng('shuffle','v5normal');


sx(1:L,1)=2*real(sigma_m_0);
sy(1:L,1)=2*imag(sigma_m_0);
sz(1:L,1)=sigma_z_0;




for m=2:M0
%     
%sample noise that comes from correlated dissipation
noise_gamma_d_diag1=is_noise*mask_l.*(randn(L,1))*sqrt(dt); % noise in diagonal basis
noise_gamma_d_diag2=is_noise*mask_l.*(randn(L,1))*sqrt(dt); % noise in diagonal basis
Gamma_l=is_noise*Gamma_l;
noiseGammaPrev_real_d=(noise_gamma_d_diag1.'*V_l).';
noiseGammaPrev_imag_d=(noise_gamma_d_diag2.'*V_l).';

% sample noise from correlated pump
noise_gamma_p_diag1=is_noise*mask_p.*(randn(L,1))*sqrt(dt); % noise in diagonal basis
noise_gamma_p_diag2=is_noise*mask_p.*(randn(L,1))*sqrt(dt); % noise in diagonal basis
Gamma_p=is_noise*Gamma_p;
noiseGammaPrev_real_p=(noise_gamma_p_diag1.'*V_p).';
noiseGammaPrev_imag_p=(noise_gamma_p_diag2.'*V_p).';
 %% previous step
    A11=-2*omega0*sy+sz.*(J_l*sy)+0.5*sz.*(Gamma_l*sx)+sz.*(J_p*sy)-0.5*sz.*(Gamma_p*sx);
    A12=2*omega0*sx-2*Omega*sz-sz.*(J_l*sx)+0.5*sz.*(Gamma_l*sy)-sz.*(J_p*sx)-0.5*sz.*(Gamma_p*sy);
    A13=2*Omega*sy+sy.*(J_l*sx)-sx.*(J_l*sy)-0.5*real(sx.*(Gamma_l*sx)+sy.*(Gamma_l*sy)) +sy.*(J_p*sx)-sx.*(J_p*sy)+0.5*real(sx.*(Gamma_p*sx)+sy.*(Gamma_p*sy));
    
     
    %% new step predictor
    sx_tilde=sx+A11*dt+noiseGammaPrev_imag_d.*sz-noiseGammaPrev_imag_p.*sz;
    sy_tilde=sy+A12*dt+noiseGammaPrev_real_d.*sz-noiseGammaPrev_real_p.*sz;
    sz_tilde=sz+A13*dt-(noiseGammaPrev_real_d.*sy+noiseGammaPrev_imag_d.*sx)+(noiseGammaPrev_real_p.*sy+noiseGammaPrev_imag_p.*sx); 

   
    for i=1:n_iter
        

    A11_tilde=-2*omega0*sy_tilde+sz_tilde.*(J_l*sy_tilde)+0.5*sz_tilde.*(Gamma_l*sx_tilde)+sz_tilde.*(J_p*sy_tilde)-0.5*sz_tilde.*(Gamma_p*sx_tilde);
    A12_tilde=2*omega0*sx_tilde-2*Omega*sz_tilde-sz_tilde.*(J_l*sx_tilde)+0.5*sz_tilde.*(Gamma_l*sy_tilde)-sz_tilde.*(J_p*sx_tilde)-0.5*sz_tilde.*(Gamma_p*sy_tilde);
    A13_tilde=2*Omega*sy_tilde+sy_tilde.*(J_l*sx_tilde)-sx_tilde.*(J_l*sy_tilde)-0.5*real(sx_tilde.*(Gamma_l*sx_tilde)+sy_tilde.*(Gamma_l*sy_tilde)) +sy_tilde.*(J_p*sx_tilde)-sx_tilde.*(J_p*sy_tilde)+0.5*real(sx_tilde.*(Gamma_p*sx_tilde)+sy_tilde.*(Gamma_p*sy_tilde));
 
    
    %% new step predictor
    sx_new=sx+0.5*(A11+A11_tilde)*dt+0.5*noiseGammaPrev_imag_d.*(sz+sz_tilde)-0.5*noiseGammaPrev_imag_p.*(sz+sz_tilde);
    sy_new=sy+0.5*(A12+A12_tilde)*dt+0.5*noiseGammaPrev_real_d.*(sz+sz_tilde)-0.5*noiseGammaPrev_real_p.*(sz+sz_tilde);
    sz_new=sz+0.5*(A13+A13_tilde)*dt-0.5*(noiseGammaPrev_real_d.*(sy+sy_tilde)+noiseGammaPrev_imag_d.*(sx+sx_tilde))+ 0.5*(noiseGammaPrev_real_p.*(sy+sy_tilde)+noiseGammaPrev_imag_p.*(sx+sx_tilde)); 

    
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

    

    sx=sx_new;
    sy=sy_new;
    sz=sz_new;

   
end
S_m= 0.5*(sx_new-1i*sy_new);
S_z= sz_new;
end