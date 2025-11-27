function [S_m,S_z]= solver_spins_corr_emission_real_loop_prep(sigma_m_0,sigma_z_0,J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p)
% this code does stochastic evolution and returns finale state at time
% tau_0 . In principle, one can also save intermediate values of S_m, S_z,
% also t as it is done in solver_spins_corr_emission_real_loop_par(). 

parameters_spin;



%% preparation
rng('shuffle','v5normal');

sx(1:N,1:ncopies)=2*real(sigma_m_0);
sy(1:N,1:ncopies)=2*imag(sigma_m_0);
sz(1:N,1:ncopies)=sigma_z_0;


mask_l=is_noise*mask_l;
mask_p=is_noise*mask_p;

ng_d1 = randn(N, ncopies, Knoise, 'like', sx) .* mask_l * sqrt(dt);
    ng_d2 = randn(N, ncopies, Knoise, 'like', sx) .* mask_l * sqrt(dt);

    ng_p1 = randn(N, ncopies, Knoise, 'like', sx) .* mask_p * sqrt(dt);
    ng_p2 = randn(N, ncopies, Knoise, 'like', sx) .* mask_p * sqrt(dt);

    % Transform noise to real space:
    % Each slice: (N×ncopies)' * (N×N) = (ncopies×N),
    % then transpose back to N×ncopies
    noise_real_d = pagemtimes(permute(ng_d1,[2 1 3]), V_l); 
    noise_real_d = permute(noise_real_d,[2 1 3]);
    noise_imag_d = pagemtimes(permute(ng_d2,[2 1 3]), V_l); 
    noise_imag_d = permute(noise_imag_d,[2 1 3]);

    noise_real_p = pagemtimes(permute(ng_p1,[2 1 3]), V_p); 
    noise_real_p = permute(noise_real_p,[2 1 3]);
    noise_imag_p = pagemtimes(permute(ng_p2,[2 1 3]), V_p); 
    noise_imag_p = permute(noise_imag_p,[2 1 3]);

knoise=1;

 Aa11=@(x,y,z)-2*omega0*y+z.*(J_l*y)+0.5*z.*(Gamma_l*x)+z.*(J_p*y)-0.5*z.*(Gamma_p*x);
    Aa12=@(x,y,z)2*omega0*x-2*Omega*z-z.*(J_l*x)+0.5*z.*(Gamma_l*y)-z.*(J_p*x)-0.5*z.*(Gamma_p*y);
    Aa13=@(x,y,z)2*Omega*y+y.*(J_l*x)-x.*(J_l*y)-0.5*real(x.*(Gamma_l*x)+y.*(Gamma_l*y)) +y.*(J_p*x)-x.*(J_p*y)+0.5*real(x.*(Gamma_p*x)+y.*(Gamma_p*y));


for m=2:M0


%     
%sample noise that comes from correlated dissipation
% if DO_PARALLEL
% noise_gamma_d_diag1=is_noise*mask_l.*(randn(N,ncopies,'gpuArray'))*sqrt(dt); % noise in diagonal basis
% noise_gamma_d_diag2=is_noise*mask_l.*(randn(N,ncopies,'gpuArray'))*sqrt(dt); % noise in diagonal basis
% Gamma_l=is_noise*Gamma_l;
% noiseGammaPrev_real_d=(noise_gamma_d_diag1.'*V_l).';
% noiseGammaPrev_imag_d=(noise_gamma_d_diag2.'*V_l).';
% 
% % sample noise from correlated pump
% noise_gamma_p_diag1=is_noise*mask_p.*(randn(N,ncopies,'gpuArray'))*sqrt(dt); % noise in diagonal basis
% noise_gamma_p_diag2=is_noise*mask_p.*(randn(N,ncopies,'gpuArray'))*sqrt(dt); % noise in diagonal basis
% Gamma_p=is_noise*Gamma_p;
% noiseGammaPrev_real_p=(noise_gamma_p_diag1.'*V_p).';
% noiseGammaPrev_imag_p=(noise_gamma_p_diag2.'*V_p).';
% else 
%     noise_gamma_d_diag1=is_noise*mask_l.*(randn(N,ncopies))*sqrt(dt); % noise in diagonal basis
% noise_gamma_d_diag2=is_noise*mask_l.*(randn(N,ncopies))*sqrt(dt); % noise in diagonal basis
% Gamma_l=is_noise*Gamma_l;
% noiseGammaPrev_real_d=(noise_gamma_d_diag1.'*V_l).';
% noiseGammaPrev_imag_d=(noise_gamma_d_diag2.'*V_l).';
% 
% % sample noise from correlated pump
% noise_gamma_p_diag1=is_noise*mask_p.*(randn(N,ncopies))*sqrt(dt); % noise in diagonal basis
% noise_gamma_p_diag2=is_noise*mask_p.*(randn(N,ncopies))*sqrt(dt); % noise in diagonal basis
% Gamma_p=is_noise*Gamma_p;
% noiseGammaPrev_real_p=(noise_gamma_p_diag1.'*V_p).';
% noiseGammaPrev_imag_p=(noise_gamma_p_diag2.'*V_p).';
% end

 noiseGammaPrev_real_d(1:N,1:ncopies)=noise_real_d(:,:,knoise);
    noiseGammaPrev_imag_d(1:N,1:ncopies)=noise_imag_d(:,:,knoise);
    noiseGammaPrev_real_p(1:N,1:ncopies)=noise_real_p(:,:,knoise);
    noiseGammaPrev_imag_p(1:N,1:ncopies)=noise_imag_p(:,:,knoise);

    knoise=knoise+1;

    if knoise==Knoise+1
 ng_d1 = randn(N, ncopies, Knoise, 'like', sx) .* mask_l * sqrt(dt);
    ng_d2 = randn(N, ncopies, Knoise, 'like', sx) .* mask_l * sqrt(dt);

    ng_p1 = randn(N, ncopies, Knoise, 'like', sx) .* mask_p * sqrt(dt);
    ng_p2 = randn(N, ncopies, Knoise, 'like', sx) .* mask_p * sqrt(dt);

    % Transform noise to real space:
    % Each slice: (N×ncopies)' * (N×N) = (ncopies×N),
    % then transpose back to N×ncopies
    noise_real_d = pagemtimes(permute(ng_d1,[2 1 3]), V_l); 
    noise_real_d = permute(noise_real_d,[2 1 3]);
    noise_imag_d = pagemtimes(permute(ng_d2,[2 1 3]), V_l); 
    noise_imag_d = permute(noise_imag_d,[2 1 3]);

    noise_real_p = pagemtimes(permute(ng_p1,[2 1 3]), V_p); 
    noise_real_p = permute(noise_real_p,[2 1 3]);
    noise_imag_p = pagemtimes(permute(ng_p2,[2 1 3]), V_p); 
    noise_imag_p = permute(noise_imag_p,[2 1 3]);

knoise=1;

    end
 %% previous step
      A11=Aa11(sx,sy,sz);
    A12=Aa12(sx,sy,sz);
    A13=Aa13(sx,sy,sz);
    % dissipative dynamics for shifted initial conditions
    % A11_p=-2*omega0*sy_p+sz_p.*(J_l*sy_p)+0.5*sz_p.*(Gamma_l*sx_p)+sz_p.*(J_p*sy_p)-0.5*sz_p.*(Gamma_p*sx_p);
    % A12_p=2*omega0*sx_p-2*Omega*sz_p-sz_p.*(J_l*sx_p)+0.5*sz_p.*(Gamma_l*sy_p)-sz_p.*(J_p*sx_p)-0.5*sz_p.*(Gamma_p*sy_p);
    % A13_p=2*Omega*sy_p+sy_p.*(J_l*sx_p)-sx_p.*(J_l*sy_p)-0.5*real(sx_p.*(Gamma_l*sx_p)+sy_p.*(Gamma_l*sy_p)) +sy_p.*(J_p*sx_p)-sx_p.*(J_p*sy_p)+0.5*real(sx_p.*(Gamma_p*sx_p)+sy_p.*(Gamma_p*sy_p));
  
    %% new step predictor
    sx_tilde=sx+A11*dt+noiseGammaPrev_imag_d.*sz-noiseGammaPrev_imag_p.*sz;
    sy_tilde=sy+A12*dt+noiseGammaPrev_real_d.*sz-noiseGammaPrev_real_p.*sz;
    sz_tilde=sz+A13*dt-(noiseGammaPrev_real_d.*sy+noiseGammaPrev_imag_d.*sx)+(noiseGammaPrev_real_p.*sy+noiseGammaPrev_imag_p.*sx); 
    % predictor for shifted initial conditions
    for i=1:n_iter
        
        sx_mid=0.5*(sx+sx_tilde);
        sy_mid=0.5*(sy+sy_tilde);
        sz_mid=0.5*(sz+sz_tilde);
% dissipative dynamics for estimated state at time m
    A11_tilde=Aa11(sx_mid,sy_mid,sz_mid);
    A12_tilde=Aa12(sx_mid,sy_mid,sz_mid);
    A13_tilde=Aa13(sx_mid,sy_mid,sz_mid);

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